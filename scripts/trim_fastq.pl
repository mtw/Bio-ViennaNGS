#!/usr/bin/env perl
# Last changed Time-stamp: <2014-12-10 13:10:50 mtw>
# AUTHOR: Joerg Fallmann <joerg.fallmann@univie.ac.at>

###############
###Use stuff
###############
use warnings;
use strict;
use PerlIO::gzip;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;

###############
###Variables
###############
my $VERBOSE = 0;
my ($file, $five, $three);

###############
###Command Line Options
###############

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
					   "fastq|f=s" => \$file,
					   "up|u=s"    => \$five,
					   "down|d=s"  => \$three,
					   "man"       => sub{pod2usage(-verbose => 2)},
					   "help|h"    => sub{pod2usage(-verbose => 1)},
					   "verbose"   => sub{ $VERBOSE++ }
					  );

###############
###MAIN
###############

unless ( $file && $five && $three) {
  warn "Please provide the name of the fastq via the -f option and definde the trim length via -u  and -d";
  pod2usage(-verbose => 0);
}

unless ( -f $file) {
  warn "Could not find input file $file given via -f option";
  pod2usage(-verbose => 0);
}

open(IN,"<:gzip(autopop)",$file)|| die "Can't open file: $!"; 
my $x=1;
while(<IN>){
	my $string="$_";
	chomp($string);
	$string =~ s/\cn/\\cn\cn/gs;
	if( (++$x)%2 ){
	    my $i = substr($string,$three,-$five);
	    print "$i\n"
	}
	else {
	    print "$string\n";
	}
}
close(IN);

__END__

=head1 NAME

trim_fastq.pl - Trims sequence and quality string of fastq files on the fly

=head1 SYNOPSIS

trim_fastq.pl [-f I<FILE>] [-u I<INTEGER>] [-d I<INTEGER>]
[options]

=head1 DESCRIPTION

This program trims the sequence and qualitystring fields from a fastq
file by user defined length, for example to allow re-mapping if
mapping quality is not sufficient without trimming.


=head1 OPTIONS

=over 

=item B<-f>

Fastq file for trimming

=item B<-u>

Number of nucleotides to trim from read/qualitystring start

=item B<-d>

Number of nucleotides to trim from read/qualitystring end

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut
