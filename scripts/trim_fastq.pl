#!/usr/bin/env perl
# Last changed Time-stamp: <2016-04-21 15:52:48 mtw>
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
# *  You should have received a copy of the GNU General Public License},
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *  This copyright notice MUST APPEAR in all copies of the script!
# ***********************************************************************

use warnings;
use strict;
use PerlIO::gzip;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($file, $five, $three, $a5, $str, $count, $last);
my $len = 0;
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
					   "fastq|f=s" => \$file,
					   "up|u=s"    => \$five,
					   "down|d=s"  => \$three,
					   "a5=s"      => \$a5,
					   "man"       => sub{pod2usage(-verbose => 2)},
					   "help|h"    => sub{pod2usage(-verbose => 1)}
					  );

unless ( $file ) {
  warn "Please provide the name of the fastq via the -f option";
  pod2usage(-verbose => 0);
}

unless ($a5) {
  unless ( $five && $three ){
    warn "Please define trim length via -u  and -d";
    pod2usage(-verbose => 0);
  }
}
else {
  if ( $five && $three ){
    warn "The -a5 and -u/-d options are mutually exclusive. Please do not use -u/-d when trimming a 5' adapter sequence";
    pod2usage(-verbose => 0);
  }
  $len = length($a5);
}

unless ( -f $file) {
  warn "Could not find input file $file given via -f option";
  pod2usage(-verbose => 0);
}

$count = 0;
$len = length($a5);

open(IN,"<:gzip(autopop)",$file)|| die "Can't open file: $!";
my $x=1;
while(<IN>){
  my $string="$_";
  chomp($string);
  $string =~ s/\cn/\\cn\cn/gs;
  if( (++$x)%2 ){
  #  print "> $x $string\n";
    if( $a5){
      my $pattern = qr/$a5/;
      if ($string =~ m/^$pattern/){
	$count++;
	$str = substr($string,$len);
	$last = $x;
      }
      if ($x == $last + 2){ # trim quality string for previously trimmed read
	$str = substr($string,$len);
	$last = -1;
      }
      print "$x $str\n";
    }
    else{
      my $i = substr($string,$three,-$five);
      print "$i\n";
    }
  }
  else {
    print "* $string\n";
  }
}
close(IN);

print "found adapter $a5 $count times\n"

__END__

=head1 NAME

trim_fastq.pl - Trim sequence and quality string of fastq files on the fly

=head1 SYNOPSIS

trim_fastq.pl [-f I<FILE>] [-u I<INTEGER>] [-d I<INTEGER>] [-a5 I<STRING>]
[options]

=head1 DESCRIPTION

This program trims sequence and qualitystring fields from a fastq
file, either by a number of nucleotides, or by removing a user defined
5' adapter sequence. The typical usage scenario is to allow re-mapping
if mapping quality is not sufficient without trimming.


=head1 OPTIONS

=over1

=item B<--fastq -f>

Fastq file for trimming

=item B<--up -u>

Number of nucleotides to trim from read / quality string start

=item B<--down -d>

Number of nucleotides to trim from read / quality string end

=item B<--a5>

Sequence of 5' adapter to trim from the left side of the reads. This
option and -u / -d are mutually exclusive

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHORS

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>
Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
