
=head1 NAME

Bio::Metabolic::Network - Perl extension for biochemical reaction networks

=head1 SYNOPSIS

  use Bio::Metabolic::Network;

  my $net = Bio::Metabolic::Network->new($reaction1, $reaction2, ... );


=head1 DESCRIPTION

This class implements objects representing biochemical reaction networks.
A reaction network is defined a number of biochemical reactions.

=head2 EXPORT

None

=head2 OVERLOADED OPERATORS

  String Conversion
    $string = "$network";
    print "\$network = '$network'\n";

  Comparison
    if ($network1 <= $network2)...


=head1 AUTHOR

Oliver EbenhÂŽöh, oliver.ebenhoeh@rz.hu-berlin.de

=head1 SEE ALSO

Bio::Metabolic Bio::Metabolic::Substrate Bio::Metabolic::Substrate::Cluster Bio::Metabolic::Reaction.

=cut

package Bio::Metabolic::Network;

require 5.005_62;
use strict;
use warnings;

require Exporter;

use Bio::Metabolic::Substrate;
use Bio::Metabolic::Substrate::Cluster;
use PDL;
use PDL::Matrix;

use Math::Symbolic;
use Math::Symbolic::VectorCalculus;

#use PDL::Matrix::Extras;

use Carp;

use overload
  "\"\"" => \&network_to_string,
  "+"    => \&add_networks,
  "<="   => \&is_in,
  ">="   => sub {
    my $n1 = shift;
    my $n2 = shift;
    return ( $n2 <= $n1 );
  },
  "==" => sub {
    my $n1 = shift;
    my $n2 = shift;
    return ( $n1 <= $n2 && $n2 <= $n1 );
  };

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Bio::Metabolic::Network ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = (
    'all' => [
        qw(

          )
    ]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(

);
our $VERSION = '0.06';

our %OutputFormat = (
    'substrate' => "%20s",
    'entry'     => "%5d"
);

# Below is stub documentation for your module. You better edit it!

=head1 METHODS

=head2 Constructor new

Returns a new Bio::Metabolic::Network object.
Passed arguments must be Bio::Metabolic::Reaction objects.

Every network object is associated with a matrx, the stoichiometric matrix. 
This matrix is defined by the reactions and gets determined upon creation.

=cut

sub new {
    my $pkg       = shift;
    my @reactions = @_ ? $_[0] =~ /ARRAY/ ? @{ $_[0] } : @_ : ();

    #  my $reactions = @_ ? shift : [ ];

    #  $reactions = [$reactions, @_] unless ref($reactions) =~ /ARRAY/;

    my $new_network = bless { reactions => \@reactions } => $pkg;
    $new_network->new_matrix;

    return $new_network;
}

=head2 Method copy

Returns a clone of the network.
However, the references to the reactions point to exactly the same reactions.
If one gets modified, it effects all networks with that reaction.

=cut

sub copy {
    my $orig = shift;
    return ref($orig)->new( $orig->reactions );
}

=head2 Method reactions

Returns an arrayref of the reactions.

=cut

sub reactions {
    return shift->{'reactions'};
}

=head2 Method has_reaction

Argument is a reaction. Returns 1 if the network contains the raction, 0 otherwise.

=cut

sub has_reaction {
    my ( $network, $reaction ) = @_;

    my $nr    = 0;
    my @rlist = @{ $network->reactions };

    #  print "rlist has ".eval(@rlist)." elements\n";
    foreach my $netr (@rlist) {
        $nr++ if ( $netr == $reaction );
    }

    return $nr;
}

=head2 Method add_reaction

Argument is a Bio::Metabolic::Reaction. 
Alters the object in-line, adding the reaction to the list.

=cut

sub add_reaction {
    my ( $network, $reaction ) = @_;

    unless ( $network->has_reaction($reaction) ) {
        push( @{ $network->reactions }, $reaction );
        $network->new_matrix;
    }
}

=head2 Method remove_reaction

Argument is a Bio::Metabolic::Reaction. 
Altering the object in-line, removeing the reaction from the network.
Return undef if network did not have the reaction, 1 on success.

=cut

sub remove_reaction {
    my ( $network, $reaction ) = @_;

    return undef unless $network->has_reaction($reaction);

    my $cut  = 0;
    my $netr = $network->reactions;
    my $cnt  = 0;
    while ( $cnt < @$netr ) {
        if ( $netr->[$cnt] == $reaction ) {
            splice( @$netr, $cnt, 1 );
            $network->new_matrix;
            $cut++;
        }
        else {
            $cnt++;
        }
    }

    #  print "remove_reaction: $cut removed.\n";
    return $cut;
}

=head2 Method network_to_string

Returns a string representation of the network

=cut

sub network_to_string {
    my $network = shift;

    my @rlist  = @{ $network->reactions };
    my $nr     = @rlist;
    my $retstr = "$nr reactions:\n";
    foreach my $r (@rlist) {
        $retstr .= $r . "-------------------------------------------\n";
    }
    return $retstr;
}

=head2 Method add_networks

Adds an arbitrary number of networks, returning a new object containing all reactions that are
contained at least present in one of the networks.

=cut

sub add_networks {
    my @nets = @_;

    # this is due to an extra value passed by the overload Module
    pop(@nets) if ( ref( $nets[ @nets - 1 ] ) ne ref( $nets[0] ) );

    croak("add_network needs at least one network!") if @nets == 0;

    my $newnet = ref( $nets[0] )->new;

    foreach my $network (@nets) {
        foreach my $reaction ( @{ $network->reactions } ) {
            $newnet->add_reaction($reaction);
        }
    }

    return $newnet;
}

=head2 Method is_in

$net1->is_in($net2) Returns 1 if all reactions in $net1 also occur in $net2,
i.e. if $net1 is a subnetwork of $net2.

=cut

sub is_in {
    my ( $net1, $net2 ) = @_;

    foreach my $reaction ( @{ $net1->reactions } ) {
        return 0 unless $net2->has_reaction($reaction);
    }

    return 1;
}

=head2 Method dist

Provides a distance measure between networks. Returns the number of reactions that are
different in the two networks.

=cut

sub dist {
    my ( $net1, $net2 ) = @_;

    my @r1 = @{ $net1->reactions };
    my @r2 = @{ $net2->reactions };

    my $dist = @r1 > @r2 ? @r1 : @r2;
    foreach my $reaction (@r1) {
        $dist-- if $net2->has_reaction($reaction);
    }

    return $dist;
}

=head2 Method substrates

Returns a Bio::Metabolic::Substrate::Cluster containing all substrates participating in at least
one reaction.
=cut

sub substrates {
    my $network = shift;

    my $cluster = Bio::Metabolic::Substrate::Cluster->new;

    return $cluster->add_clusters(
        map( $_->substrates, @{ $network->reactions } ) );
}

=head2 Method matrix

Returns the stoichiometric matrix of the network as a PD::Matrix object.

=cut

sub matrix {
    my $network = shift;

    $network->{matrix} = shift if @_;
    return $network->{matrix};
}

=head2 Method new_matrix

determines the stoichiometric matrix of the network defined by its reactions.

=cut

sub new_matrix {
    my $network = shift;

    my $reactions  = $network->reactions;
    my $substrates = $network->substrates;

    my @slist = $substrates->list;
    my $cols  = @$reactions;
    my $rows  = @slist;

    #  croak("cannot create matrix from nothing!") if $cols == 0 || $rows == 0;
    if ( $cols == 0 || $rows == 0 ) {

        # changed 7.11.02
        #    $network->matrix(PDL::Mat->new(1,1));
        $network->matrix( PDL::Matrix->null );
        return undef;
    }

    # changed 7.11.02
    #  my $matrix = PDL::Mat->new($cols,$rows);
    my $matrix = mzeroes( $rows, $cols );

    for ( my $r = 0 ; $r < $cols ; $r++ ) {
        foreach my $dir ( -1, 1 ) {
            foreach my $substrate ( $reactions->[$r]->dir($dir)->list ) {
                my $s = $substrates->which($substrate);

             #	print "setting ($r,$s) to ".eval($matrix->at($r,$s) + $dir)."\n";
             #	set $matrix, $r, $s, $matrix->at($r,$s) + $dir;

                # changed 7.11.02
                #	$matrix->set($r,$s,$matrix->at($r,$s) + $dir);
                $matrix->set( $s, $r, $matrix->at( $s, $r ) + $dir );
            }
        }
    }

    #  print "matrix:\n$matrix";

    #  return bless {'matrix' => $matrix,
    #		'substrates' => $substrates,
    #		'reactions' => $reactions
    #	       }, $pkg;
    $network->matrix($matrix);
}

=head2 Method print_matrix

Prints the matrix in a way that describes which substrates are associated with which rows
and which reactions with which columns.

=cut

sub print_matrix {
    my $network = shift;

    my $m     = $network->matrix;
    my @slist = $network->substrates->list;

    # changed 7.11.02
    #  my ($cols,$rows) = $m->dims;
    my ( $rows, $cols ) = $m->mdims;

    my $retstr = sprintf( $OutputFormat{'substrate'}, "" );
    for ( my $r = 0 ; $r < $cols ; $r++ ) {
        $retstr .= sprintf( $OutputFormat{'entry'}, $r );
    }
    $retstr .= "\n";

    for ( my $s = 0 ; $s < $rows ; $s++ ) {
        $retstr .= sprintf( $OutputFormat{'substrate'}, $slist[$s] . ": [" );
        for ( my $r = 0 ; $r < $cols ; $r++ ) {

            # changed 7.11.02
            #      $retstr .= sprintf($OutputFormat{'entry'},$m->at($r,$s));
            $retstr .= sprintf( $OutputFormat{'entry'}, $m->at( $s, $r ) );
        }
        $retstr .= "]\n";
    }

    return $retstr;
}

=begin comment

_time_derivative_by_substrate_number returns the temporal change of a substrate
(specified by an index). This change is determined by teh reaction rates.

=end comment

=cut

sub _time_derivative_by_substrate_number {
    my $network = shift;
    my $s       = shift;

    my $matrix    = $network->matrix;
    my $reactions = $network->reactions;

    my $sterm = Math::Symbolic::Constant->zero();

    my @rcols = PDL::which( $matrix->slice("($s),:") != 0 )->list;
    foreach my $r (@rcols) {
        my $rterm = Math::Symbolic::Operator->new(
            '*',
            Math::Symbolic::Constant->new( $matrix->at( $s, $r ) ),
            $reactions->[$r]->rate
        );

#    my $rterm = Symbolic::Function::Operator::Mult->new([Symbolic::Function::Constant->new($matrix->at($s,$r)),
#							 $reactions->[$r]->rate])->simplify;
        $sterm =
          Math::Symbolic::Operator->new( '+', $rterm->simplify, $sterm )
          ->simplify;

#    $sterm = Symbolic::Function::Operator::Add->new([$rterm,$sterm])->simplify;
    }

    return $sterm;
}

=head2 Method time_derivative

This method returns a Math::Symbolic tree representing a function which determines the temporal
change of the concentration of the substrate passed as the argument.
This function is determined by the reaction rates.

=cut

sub time_derivative {
    my $network   = shift;
    my $substrate = shift;

    my $s = $network->substrates->which($substrate);

    return
      defined $s ? $network->_time_derivative_by_substrate_number($s) : undef;
}

=head2 Method ODEs

This method returns an array of Math::Symbolic trees. One for each substrate.
The optional argument can be a Bio::Metabolic::Substrate::Cluster, an arrayref or an array
of Bio::Metabolic::Substrate objects.
If no argument is specified, it defaults to what the method substrates() returns.

=cut

sub ODEs {
    my $network = shift;

    my @sindices   = ();
    my @substrates = ();
    if (@_) {
        if ( ref( $_[0] ) eq 'Bio::Metabolic::Substrate::Cluster' ) {
            @substrates = shift->list;
        }
        elsif ( ref( $_[0] ) eq 'ARRAY' ) {
            @substrates = @{ shift() };
        }
        else {
            @substrates = @_;
        }
    }
    else {
        @substrates = $network->substrates->list;
    }

    #  if (@_) {
    #    my @substrates = shift->list;
    while ( my $sub = shift(@substrates) ) {
        push( @sindices, $network->substrates->which($sub) );
    }

    my @time_derivatives = ();

    for ( my $s = 0 ; $s < @sindices ; $s++ ) {
        push( @time_derivatives,
            $network->_time_derivative_by_substrate_number( $sindices[$s] ) );

        #    push (@time_derivatives,$sterm);
    }

    return @time_derivatives;
}

=head2 Method mfile

Dumps a strong which can be used as an mfile for Matlab.
The optional argument can be a Bio::Metabolic::Substrate::Cluster, an arrayref or an array
of Bio::Metabolic::Substrate objects.
If no argument is specified, it defaults to what the method substrates() returns.

The corresponding substrate concentrations are considered to be integration variables.

All parameters within the participating reactions must have a value.
The method croaks if one parameter value is not defined.

=cut

sub mfile {
    my $network = shift;

    #  my $substrates = shift;

    my @substrates = ();
    if (@_) {
        if ( ref( $_[0] ) eq 'Bio::Metabolic::Substrate::Cluster' ) {
            @substrates = shift->list;
        }
        elsif ( ref( $_[0] ) eq 'ARRAY' ) {
            @substrates = @{ shift() };
        }
        else {
            @substrates = @_;
        }
    }
    else {
        @substrates = $network->substrates->list;
    }

#  my @substrates = ref($_[0]) ? ref($_[0]) eq 'Bio::Metabolic::Substrate::Cluster' ? $_[0]->list : @{$_[0]} : @_;

    my @odes = $network->ODEs(@substrates);

    #  my @odes = &$odesub($network, @substrates);

    #  my $varlist = Math::Symbolic::VectorCalculus::_combined_signature(@odes);

    my %varchecklist = map ( ( $_, 1 ),
        @{ Math::Symbolic::VectorCalculus::_combined_signature(@odes) } );

    my %parameters = ();    # collect all parameter values;
    foreach my $reaction ( @{ $network->reactions } ) {
        foreach my $param ( keys( %{ $reaction->parameters } ) ) {
            $parameters{ $reaction->parameter($param)->{name} } =
              $reaction->parameter($param)->value();
        }
    }

    my $mfile = "function f=func(t,y)\n";

    my @varnames = map ( $_->name, @substrates );
    for ( my $vnr = 0 ; $vnr < @varnames ; $vnr++ ) {
        $mfile .= $varnames[$vnr] . "=y(" . eval( $vnr + 1 ) . ");\n";
        delete( $varchecklist{ $varnames[$vnr] } )
          if defined( $varchecklist{ $varnames[$vnr] } );
    }

    foreach my $param ( keys(%varchecklist) ) {
        croak("undefined parameter $param in method mfile")
          unless defined $parameters{$param};
        $mfile .= $param . "=" . $parameters{$param} . "\n";

        foreach my $ode (@odes) {
            $ode->implement( $param => $parameters{$param} );
        }

        delete( $varchecklist{$param} ) if defined( $varchecklist{$param} );
    }

    for ( my $snr = 0 ; $snr < @odes ; $snr++ ) {
        my ( $codeline, $leftovers ) =
          Math::Symbolic::Compiler->compile_to_code( $odes[$snr], \@varnames );
        croak("cannot handle leftover trees in method mfile") if @$leftovers;

        for ( my $vnr = 0 ; $vnr < @varnames ; $vnr++ ) {
            $codeline =~ s/\$_\[$vnr\]/$varnames[$vnr]/g;
        }

        $mfile .= "f(" . eval( $snr + 1 ) . ",1)=" . $codeline . ";\n";
    }

    return $mfile;
}


1;
__END__
# construction site!!!!!!
sub can_convert {
  my $network = shift;

    my @substrates = ();
    if (@_) {
        if ( ref( $_[0] ) eq 'Bio::Metabolic::Substrate::Cluster' ) {
            @substrates = shift->list;
        }
        elsif ( ref( $_[0] ) eq 'ARRAY' ) {
            @substrates = @{ shift() };
        }
        else {
            @substrates = @_;
        }
    }

  my $slist = $network->substrates;
  foreach my $ext (@substrates) {
      push (@ex_indices, $slist->which($ext)) if defined $slist->which($ext);
  }


#  print "deleting rows (".join(',',@ex_indices).")\n";
  my $reduced = $network->matrix->copy;
  $reduced->delrows(@ex_indices);
#  print "reduced matrix is : $reduced\n";

  my $kernel = $reduced->kernel();
#  print "kernel is $kernel\n";
  return undef if $kernel->isempty;
#  my ($kerneldim,$nor) = $kernel->dims();

  my $convert = $network->matrix x $kernel;

#  print "convert: $convert\n";
  foreach $ext (@ex_indices) {
# changed 8.11.02
#    my $res = $convert->slice(":,($ext)");
    my $res = $convert->slice("($ext),:");
    return undef if ($res->where($res != 0)->isempty);
  }

  return $kernel;
}

sub is_elementary {
    my $net = shift;

    my $kernel = $net->can_convert(@_);

    return undef unless defined($kernel);

    my @kdims = $kernel->mdims;

    return undef if $kdims[1] != 1;

    return undef if which($kernel->slice(":,(0)")==0)->nelem > 0;

    my $conversion = $net->matrix x $kernel;
    my @cdims = $conversion->mdims;

    return undef if which($conversion->slice(":,(0)")==0)->nelem == $cdims[0];

    return $kernel;
}

#sub conservation_rules {
#  my $network = shift;

#  my $externals = ref($_[0]) =~ /Bio::Metabolic::Substrate::Cluster/ ? $_[0] :
#    Bio::Metabolic::Substrate::Cluster->new(@_);

#  my $st = $network->matrix->copy;
#  my $slist = $network->substrates;

#  my @ext_indeces = ();
#  foreach my $ext ($externals->list) {
#    push (@ext_indeces,$slist->which($ext));
#  }

#  $slist->remove_substrates($externals->list);

#  my $redt = $st->cutrows(@ext_indeces)->xchg(0,1);

#  my $crules = $redt->kernel;
##  print "kernel:\n $crules\n";
#  my @crdims = $crules->dims;
##  print "dimensions: (".join(',',@crdims).")\n";

#  my @rules = ();
#  for (my $i=0;$i<$crdims[0];$i++) {
##    print "generating rule with substrates: ".join(',',$slist->list)."\n";
#    push (@rules, Bio::Metabolic::ConservationRule->get_or_new($slist,$crules->slice("($i),:")));
#  }

#  return wantarray ? @rules : \@rules;
#}

							
sub fix_influx {
  my $network = shift;
  my $substrate = shift;
  my $influx = shift;

#  print "infinte source is: $INFINITE_SOURCE\n";
#  print "ref: ".ref($INFINITE_SOURCE)."\n";
#  print "substrate is: $substrate\n";
#  print "ref: ".ref($substrate)."\n";

  my $in_reaction = Bio::Metabolic::Reaction->get_or_new("influx_".$substrate->name,[$INFINITE_SOURCE],[$substrate]);

#  print "reaction created: $in_reaction\n";

  $in_reaction->constant(-1)->set($influx);
  $in_reaction->constant(1)->set(0);

  $network->add_reaction($in_reaction);
}

sub fix_outflux {
  my $network = shift;
  my $substrate = shift;
  my $outflux = shift;

  my $out_reaction = Bio::Metabolic::Reaction->get_or_new("outflux_".$substrate->name,[$substrate],[$INFINITE_SOURCE]);

  $out_reaction->constant(-1)->set($outflux);
  $out_reaction->constant(1)->set(0);

  $network->add_reaction($out_reaction);
}

sub assemble {
  # call: Network->assemble([{'c'=>6},{'c'=1}], $net_db)
  # or: $net->assemble([{'c'=>6},{'c'=1}], $net_db)
  my $network = shift;
  my $externals = shift;
  my $database = shift;

  if (!ref($network)) {
    my $pkg = $network;
    $network = $pkg->new();
    my @rlist = $database->fetch_all_reactions;
    $network->add_reaction($rlist[int(rand(@rlist))]);
  }

  return $network if defined($network->can_convert($externals));

  my $substrates = $network->substrates;

  my ($cols,$rows) = $network->matrix->dims;
  my @posext = ();
  for (my $r=0;$r<$rows;$r++) {
    if (which($network->matrix->slice(":,($r)") != 0)->nelem != 0) {
      push (@posext, $substrates->[$r]);
    }
  }

  my @rchoice;
  if (@posext == 0) {
    @rchoice = $database->fetch_all_reactions;
  } else {
    my $rsub = $posext[int(rand(@posext))];
    @rchoice = $database->fetch_reaction_with_substrate($rsub);
  }

  $network->add_reaction($rchoice[int(rand(@rchoice))]);

  return $network->assemble($externals, $database);
}



# experimental area for drawing the motherfuckers!

  
#sub makegraphs {
#  my $network = shift;

#  my ($in_sub, $out_sub) = ref($_[0]) =~ /ARRAY/ ? @{$_[0]} : @_;


#  my $netgraphs = Bio::Metabolic::Network::Graph->new_from_network($network, $in_sub, $out_sub);
#  print "new graph generated (".eval(@$netgraphs)." designs)\n";

#  my @graphs = ();

#  $|=1;
#  foreach my $ng (@$netgraphs) {
#    print ".";
#    push (@graphs, $ng->to_GD);
#  }
#  print "\n";
#  $|=0;

#  return wantarray ? @graphs : \@graphs;
#}
