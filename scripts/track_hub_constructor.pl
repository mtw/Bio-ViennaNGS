#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-19 23:16:22 mtw>
#
# Construct UCSC genome browser Track Hub and display various
# genomic sequence annotation data within the Hub
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael T. Wolfinger <michael@wolfinger.eu>
# *                 Florian Eggenhofer <florian.eggenhofer@univie.ac.at>
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

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Path::Class;
use Bio::ViennaNGS::UCSC qw( make_assembly_hub make_track_hub );

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($track_hub_return_value,$lf);
my $logname = "log.txt";
my $genome_identifier = '-';
my $folder_in = '-';
my $dest = '.';
my $base_URL = '-';
my $chrom_size_file="-";
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("gi|g=s"       => \$genome_identifier,
					   "infolder|i=s" => \$folder_in,
					   "out|o=s"      => \$dest,
					   "baseurl|b=s"  => \$base_URL,
                                           "chromsize|c=s"  => \$chrom_size_file,
					   "man"          => sub{pod2usage(-verbose => 2)},
					   "help|h"       => sub{pod2usage(1)}
					  );

unless ($genome_identifier){
  warn "Please provide genome identifier via --g option";
  pod2usage(-verbose => 0);
}
unless (-d $folder_in){
  warn "An input folder path is required and must be gi ven via --infolder option";
  pod2usage(-verbose => 0);
}
unless ($base_URL =~ /^http/) {
  warn "Base URL must be given in full format, eg http://foo.bar.com/Hubs/";
  pod2usage(-verbose => 0);
}
unless ($dest =~ /\/$/){$dest.= "/";}
unless (-d $dest){
  mkdir $dest or die $!;
}
$lf = file($dest,$logname);

$track_hub_return_value = make_track_hub($genome_identifier,$folder_in,$dest,$base_URL,$chrom_size_file,$lf);


__END__


=head1 NAME

track_hub_constructor.pl - Build UCSC genome browser Track Hubs from
genomic sequence and annotation

=head1 SYNOPSIS

track_hub_constructor.pl [--gi I<ID>] [--infolder I<PATH>] [--out
I<PATH>] [--baseurl -I<URL>] [options]

=head1 DESCRIPTION

The UCSC genome browser offers the possibility to visualize additional
tracks for organisms that are included in the standard UCSC browser
bundle via so called 'Track Hubs'. This script constructs Track Hubs
from annotation data.

=head1 OPTIONS

=over

=item B<--gi -g>

Genome id as used in UCSC assembly hub. Must be correct, otherwise the
annotation cannot be mapped on the genome.

=item B<--infolder -i>

Directory which contains all track files in BED/bigBed format. The
resulting Track Hub will contain these files in their respective
bigFile version.

=item B<--out -o>

Destination folder for the output Track Hub.

=item  B<--baseurl -b>

BaseURL used within the Track Hub. This URL will be included verbatim
in the resulting Track Hub. It is crucial that this URl is valid, else
the resulting Track Hub will be broken.

=item B<--help -h>

Print short help.

=item B<--man>

Prints the manual page and exits.

=back

=head1 SEE ALSO

The L<UCSC Genome
Wiki|http://genomewiki.ucsc.edu/index.php/Track_Hubs> has extensive
documentation for Assembly and Track Hubs.

=head1 AUTHORS

=over

=item Florian Eggenhofer E<lt>florian.eggenhofer@univie.ac.atE<gt>

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=back
