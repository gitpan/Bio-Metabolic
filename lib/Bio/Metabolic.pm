package Bio::Metabolic;

require 5.005_62;
use strict;
use warnings;

require Exporter;

use Bio::Metabolic::Substrate;
use Bio::Metabolic::Substrate::Cluster;
use Bio::Metabolic::Reaction;
use Bio::Metabolic::Network;

#use Bio::Metabolic::Network::Graph;
#use Bio::Metabolic::NetworkDB;
#use Bio::Metabolic::ConservationRule;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Bio::Metabolic ':all';
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

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

Bio::Metabolic - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Bio::Metabolic;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Bio::Metabolic, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.


=head1 AUTHOR

A. U. Thor, a.u.thor@a.galaxy.far.far.away

=head1 SEE ALSO

perl(1).

=cut
