#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-01-16 15:34:31 mtw>
#
# Converts bed6 to bed12
#
# usage: bed6_to_bed12.pl -i <STRING>
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2016 Michael T. Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# *  This program is free software: you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *  This copyright notice MUST APPEAR in all copies of the script!
# ***********************************************************************
#

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Bio::ViennaNGS::FeatureIO;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($fio);
my $infile  = undef;
my $name   = "foobar";


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("input|i=s"   => \$infile,
					   "man"         => sub{pod2usage(-verbose => 2)},
					   "help|h"      => sub{pod2usage(1)}
					  );

unless (-f $infile){
  warn "Could not find input file provided via --input|-i option";
  pod2usage(-verbose => 0);
}

#unless (defined $pattern) {
#  warn "No search pattern/substring provided";
#  pod2usage(-verbose => 0);
#}

$fio = Bio::ViennaNGS::FeatureIO->new(
				      file => $infile,
				      filetype => 'Bed6',
				      instanceOf => 'FeatureChain',
				     );

#print Dumper($fio);
foreach my $d (@{$fio->data}){
#  print " >> processing a ".ref($d)." object\n";
  print $d->as_bed12_line(undef,undef,undef) ."\n";
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#




__END__

=head1 NAME

bed6_to_bed12.pl - Convert BED6 to BED12.

=head1 SYNOPSIS

bed6_to_bed12.pl [--input|-i I<FILE>] [options]

=head1 DESCRIPTION

This script is a simple converter for BED6 to BED12. All it does is
compute values for the additional 6 fields according to the L<BED12
specification|https://genome.ucsc.edu/FAQ/FAQformat#format1>.

=head1 OPTIONS

=over

=item B<input|i>

A BED6 file

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
