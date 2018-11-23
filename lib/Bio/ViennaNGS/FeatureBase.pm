# -*-CPerl-*-
# Last changed Time-stamp: <2018-11-23 18:01:42 mtw>

=head1 NAME

Bio::ViennaNGS::FeatureBase - A Moose Role for BED compliance

=head1 SYNOPSIS

  package MyClass;
  use Moose;

  with 'Bio::ViennaNGS::FeatureBase';


=head1 DESCRIPTION

L<Bio::ViennaNGS::FeatureBase> is a simple L<Moose::Role> which
indicates whether or not the B<start> attribute of any of the
L<Bio::ViennaNGS> Feature* classes is zero-based or one-based
(i.e. whether or not the Feature adheres to standard BED notation).

=cut

package Bio::ViennaNGS::FeatureBase;

use Moose::Role;
use Bio::ViennaNGS::Subtypes;

has 'base' => (
	       is => 'rw',
	       isa => 'Bio::ViennaNGS::ZeroOrOne',
	       default => 0,
	      );

=head1 DEPENDENCIES

=over

=item L<Moose::Role>

=back

=head1 BUGS

Please report any bugs or feature requests through the web interface
at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-ViennaNGS>. I
will be notified, and then you'll automatically be notified of
progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::ViennaNGS::FeatureBase


You can also look for information at:

=over 2

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-ViennaNGS>

=item * Search metaCPAN

L<https://metacpan.org/release/Bio-ViennaNGS>

=back

=head1 AUTHORS

Michael T. Wolfinger, C<< <michael at wolfinger.eu> >>

=head1 LICENSE AND COPYRIGHT

Copyright 2018 Michael T. Wolfinger <michael@wolfinger.eu> and
<michael.wolfinger@univie.ac.at>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

=cut


