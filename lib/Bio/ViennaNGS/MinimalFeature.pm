# -*-CPerl-*-
# Last changed Time-stamp: <2016-11-02 16:07:51 mtw>

package Bio::ViennaNGS::MinimalFeature;

use version; our $VERSION = qv('0.16');
use namespace::autoclean;

use Moose::Util::TypeConstraints;
subtype 'PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- nor ."};
no Moose::Util::TypeConstraints;

use Moose;
with 'MooseX::Clone';

extends 'Bio::ViennaNGS::FeatureInterval';

has 'strand' => (
		 is      => 'rw',
		 isa     => 'PlusOrMinus',
		 default => '.',
		 predicate => 'has_strand',
		);

no Moose;

1;

__END__


=head1 NAME

Bio::ViennaNGS::MinimalFeature - A Moose wrapper for stranded genomic
intervals.

=head1 SYNOPSIS

  use Bio::ViennaNGS::MinimalFeature;

  my $feat = Bio::ViennaNGS::MinimalFeature->new(chromosome => "chr1",
                                                 start => "1200",
                                                 end => "4300",
                                                 strand => "+",
                                                );
=head1 DESCRIPTION

This module provides an object-oriented interface for storing
elementary stranded genomic intervals characterized via chromosome,
start position, end position and strand. As such, it can be regarded a
simple wrapper for BED4 elements.

This class inherits from L<Bio::ViennaNGS::FeatureInterval> and is the
base classs for L<Bio::ViennaNGS::Feature>.

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::FeatureInterval>

=item L<Bio::ViennaNGS::Feature>

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
