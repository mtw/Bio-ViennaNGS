#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-27 14:21:53 mtw>
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
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Path::Class;
use Bio::ViennaNGS qw(uniquify_bam);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $bam_in;
my $logfile = "bam_uniq.log";
my $outdir = "./";
my @result = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("bam=s"   => \$bam_in,
					   "o=s"     => \$outdir,
					   "man"     => sub{pod2usage(-verbose => 2)},
					   "help|h"  => sub{pod2usage(1)}
					  );


unless ($bam_in =~ /^\//) {$bam_in = "./".$bam_in;}
unless (-f $bam_in){
  warn "Could not find input file $bam_in given via -bam option";
  pod2usage(-verbose => 0);
}

#TODO check if we are allowed to write to $outdir
unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}

# Make BED12 line from each splice junction
@result = uniquify_bam($bam_in,$outdir,$logfile);

__END__


=head1 NAME

bam_uniq.pl - Deconvolute BAM files into unique and multi mappers

=head1 SYNOPSIS

bam_uniq.pl [-bam I<FILE>] [-o I<DIR>]
[options]

=head1 DESCRIPTION

This program extracts uniquely and multiply mapping reads from BAM
files acoording to the NH:i SAM attribute and creates separate BAM
files for unique and multi mappers, respectively.

=head1 OPTIONS

=over

=item B<--bam>

BAM file to extract unique and multi mappers from

=item B<-o>

Output path

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

