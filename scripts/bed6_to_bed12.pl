#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2016-10-06 16:12:42 mtw>
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

motiffinda.pl - Find sequence motifs in annotated features.

=head1 SYNOPSIS

motiffinda.pl [--motif I<REGEX>] [--gff I<FILE>] [--fa I<FILE>]
[options]

=head1 DESCRIPTION

This script extracts sequence motifs from gene annotation data. The
motif can be provided as regualr expression. Optionally, only motifs
I<in frame> with the annoation are reported.

The tool returns B<all motifs> matching the search criteria as a
(mutli-)Fasta file to STDOUT. This means if teh motif is found more
than once within an annotated feature, all matches will be reported.

=head1 OPTIONS

=over

=item B<motif|m>

The motif to search for as regular expression. For technical reasons,
the regular expression must be enclosed in brackets ().

=item B<gff|g>

Genome annotation in GFF3 format

=item B<fa|f>

Reference genome in Fasta format

=item B<gbkey>

Motifs are only searched for in this feature type, aka Genbank key
(e.g. CDS or gene). Currently only one Genbank key provided via this
option will be procesed.

=item B<offset|o>

Offset for determination of frame

=item <inframe|i>

Only report motifs in current ORF


=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
