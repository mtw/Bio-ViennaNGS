# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-22 17:01:47 mtw>

package ViennaNGS::UCSC;

use Exporter;
use version; our $VERSION = qv('0.01');
use strict;
use warnings;
use Template;

our @ISA = qw(Exporter);

our @EXPORT = qw( make_assembly_hub  );

our @EXPORT_OK = ();




1;
__END__

=head1 NAME

ViennaNGS::UCSC - Perl extension for easy UCSC Genome Browser
integration.

=head1 SYNOPSIS

  use ViennaNGS::UCSC;

=head1 DESCRIPTION

ViennaNGS::UCSC is a Perl extension for managing routine tasks with the
UCSC Genome Browser. It comes with a set of utilities that serve as
reference implementation of the routines implemented in the library. All
utilities are located in the 'scripts' folder.

=head2 EXPORT

Routines: 
  make_assembly_hub

Variables:
  none


=head3 make_assembly_hub()

Documentation for this routine goes here.

=head1 SEE ALSO

perldoc ViennaNGS
perldoc ViennaNGS::AnnoC

=head1 AUTHORS

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>
Florian Eggenhofer, E<lt>egg@tbi.univie.ac.atE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


=cut
