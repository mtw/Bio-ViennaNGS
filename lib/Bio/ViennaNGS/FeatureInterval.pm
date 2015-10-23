# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-23 16:07:24 mtw>

package Bio::ViennaNGS::FeatureInterval;

use version; our $VERSION = qv('0.16_01');
use namespace::autoclean;

use Moose;

has 'chromosome' => (
		     is  => 'rw',
		     isa => 'Str',
		     required => 1,
		     predicate => 'has_chromosome',
	     );

has 'start' => (
		is      => 'rw',
		isa     => 'Int',
		required => 1,
		predicate => 'has_start',
	       );

has 'end' => (
	      is      => 'rw',
	      isa     => 'Int',
	      required => 1,
	      predicate => 'has_end',
	     );

no Moose;

1;

__END__

=head1 NAME

Bio::ViennaNGS::FeatureInterval - A Moose class for storing elementary
genomic intervals.

=head1 SYNOPSIS

  use Bio::ViennaNGS::FeatureInterval;

  my $obj = Bio::ViennaNGS::FeatureInterval->new(chromosome => "chr1",
                                                 start => "1200",
                                                 end => "4300",
                                                );

=head1 DESCRIPTION

This module provides an object-oriented interface for storing
elementary genomic intervals characterized via chromosome, start and
end position. As such, it can be regarded a simple wrapper for BED3
elements.

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
