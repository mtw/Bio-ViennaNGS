# -*-CPerl-*-
# Last changed Time-stamp: <2018-07-03 12:59:06 mtw>

package Bio::ViennaNGS::FeatureIntervalN;

use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");
use namespace::autoclean;
use Moose;

extends 'Bio::ViennaNGS::FeatureInterval';

has 'name' => (
               is  => 'rw',
               isa => 'Str',
               required => '1',
               predicate => 'has_name',
              );

no Moose;

1;

__END__


=head1 NAME

Bio::ViennaNGS::FeatureIntervalN - A Moose wrapper for named genomic
intervals.

=head1 SYNOPSIS

  use Bio::ViennaNGS::FeatureIntervalN;

  my $feat = Bio::ViennaNGS::FeatureIntervalN->new(chromosome => "chr1",
                                                   start => "1200",
                                                   end => "4300",
                                                   name => "myFeature",
                                                   );

=head1 DESCRIPTION

This module provides an object-oriented interface for storing
elementary names genomic intervals characterized via chromosome,
start position, end position and name.

This class inherits from L<Bio::ViennaNGS::FeatureInterval>.

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::FeatureInterval>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2017 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

=cut
