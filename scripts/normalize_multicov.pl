#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-14 23:42:49 mtw>
#
# Compute normalized expression data in TPM from (raw) read
# counts in multicov.
# TPM reference: Wagner et al, Theory Biosci. 131(4), pp 281-85 (2012)
#
# usage: normalize_multicov.pl -i file.multicov -rl 50 -o /output/path
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael T. Wolfinger <michael@wolfinger.eu>
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
use Data::Dumper;
use File::Basename;
use Pod::Usage;
use Cwd;
use Path::Class;
use Bio::ViennaNGS::Util qw( parse_multicov write_multicov featCount_data computeTPM);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($infile,$FC,$FC_sample,$conds,$totreads,$i,$sample,$basename,$dir,$ext,$cwd);
my $readlength  = 100;
my $dest = getcwd();;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("i=s"             => \$infile,
					   "readlength|r=s"  => \$readlength,
					   "out|o=s"         => \$dest,
					   "man"             => sub{pod2usage(-verbose => 2)},
					   "help|h"          => sub{pod2usage(1)}
					  );

unless(-f $infile){
  warn "Could not find input multicov file provided via -i option";
  pod2usage(-verbose => 0);
  }
$cwd = getcwd();
unless ($infile =~ /^\// || $infile =~ /^\.\//){$infile = file($cwd,$infile);}
unless ($dest =~ /\/$/){$dest .= "/";}
unless (-d $dest){mkdir $dest;}
($basename,$dir,$ext) = fileparse($infile,qr/\.[^.]*/);

# parse multicov file; populate @featCount
$conds = parse_multicov($infile);
$FC = featCount_data();
#print Dumper($FC);
#print "parsing $infile: $conds conditions found\n";

# extract hash for a specific sample from $FC
for ($i=0;$i<$conds;$i++){
  my $meanTPM;
  $FC_sample = @$FC[$i];
  $meanTPM = computeTPM($FC_sample, $readlength);
}

# write multicov file based on TPM data in @featCount
write_multicov("TPM", $dest, $basename);

__END__


=head1 NAME

normalize_multicov.pl - Compute normalized expression data from read
counts

=head1 SYNOPSIS

normalize_multicov.pl [-i I<FILE>] [--readlength I<INT>] [options]

=head1 DESCRIPTION

This program computes normalized expression values in Transcript per
Million (TPM) from read counts. The latter must be provided in the
format produced by the 'bedtools multicov' utility, i.e. an extended
BED6 file where each colum past the 6th lists the read counts for one
sample/condition.

=head1 OPTION

=over

=item B<-i>

Input file in 'bedtools multicov' output format, i.e. an extended BED6
file where each colum past the 6th lists the read counts for one
sample/condition.

=item B<--readlength -r>

Read length of the RNA-seq experiment.

=item B<--out -o>

Output folder.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 SEE ALSO

TPM reference: Wagner et al, Theory Biosci. 131(4), pp 281-85 (2012)

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut




