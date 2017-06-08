# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 18:50:38 michl>

package Bio::ViennaNGS::ExtFeature;

use version; our $VERSION = qv('0.17');

use Moose;
extends 'Bio::ViennaNGS::Feature';

has 'extension' => (
		    is      => 'rw',
		    isa     => 'Str',
		    default => '',
		    predicate => 'has_extension',
		   );

no Moose;

1;

__END__


=head1 NAME

Bio::ViennaNGS::ExtFeature - A Moose wrapper for extended BED6-type
genomic intervals.

=head1 SYNOPSIS

  use Bio::ViennaNGS::ExtFeature;

  my $expression = Bio::ViennaNGS::Feature->new(chromosome => "chr1",
                                                start => "1200",
                                                end => "4300",
                                                strand => "+",
                                                name => "foobar",
                                                score => "100",
                                                extension => "34:x:36 100:200",
                                               );

=head1 DESCRIPTION

This module provides an object-oriented interface for storing extended
genomic intervals characterized via chromosome, start position, end
position, name, score, strand as well as an arbitrary string. As such,
it can be regarded a simple wrapper for extended BED6 elements.

This class inherits from L<Bio::ViennaNGS::Feature>.

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::Feature>

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
