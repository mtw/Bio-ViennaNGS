# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-23 16:41:47 mtw>

package Bio::ViennaNGS::Feature;

use version; our $VERSION = qv('0.16_01');

use Moose;
with 'MooseX::Clone';
use MooseX::InstanceTracking;

extends 'Bio::ViennaNGS::MinimalFeature';

has 'name' => (
	       is  => 'rw',
	       isa => 'Str',
	       required => '1',
	       predicate => 'has_name',
	      );

has 'score' => (
		is => 'rw',
		isa => 'Value',
		default => '0',
		predicate => 'has_score',
	       );
}

no Moose;

1;

__END__


=head1 NAME

Bio::ViennaNGS::Feature - A Moose wrapper for BED6-type genomic intervals.

=head1 SYNOPSIS

  use Bio::ViennaNGS::Feature;

  my $expression = Bio::ViennaNGS::Feature->new(chromosome => "chr1",
                                                start => "1200",
                                                end => "4300",
                                                strand => "+",
                                                name => "foobar",
                                                score => "100",
                                               );

=head1 DESCRIPTION

This module provides an object-oriented interface for storing genomic
intervals characterized via chromosome, start position, end position,
name, score and strand. As such, it can be regarded a simple wrapper
for BED6 elements.

This class inherits from L<Bio::ViennaNGS::MinimalFeature> and is the
base classs for L<Bio::ViennaNGS::FeatureChain>.

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::MinimalFeature>

=item L<Bio::ViennaNGS::FeatureChain>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2015 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

=cut
