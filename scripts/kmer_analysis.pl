#!/usr/bin/env perl
# Last changed Time-stamp: <2014-12-10 13:05:43 mtw>
# AUTHOR: Joerg Fallmann <joerg.fallmann@univie.ac.at>

###############
###Use stuff
###############
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use PerlIO::gzip;
use Cwd;
use List::Util qw(sum);
use File::Path qw(make_path remove_tree);
use Bio::ViennaNGS;

###############
###Variables
###############

my $VERBOSE = 0;
my ( $dir, $odir, $file, $klength, $type );

###############
###Command Line Options
###############

pod2usage(-verbose => 0)
	unless GetOptions(
	"dir|d=s"	=> \$dir,
	"odir|o=s"	=> \$odir,
	"file|f=s"	=> \$file,
	"kmer|k=s"      => \$klength,
	"type|t=s"      => \$type,
	"help|h"	=> sub{pod2usage(-verbose => 1)},
	"man|m"		=> sub{pod2usage(-verbose => 2)},
	"verbose"	=> sub{ $VERBOSE++ }
	);

$dir = cwd() unless ($dir);
$odir = "$dir"."/Kmer_Analysis" unless $odir;
$dir =~ s/ //g;
$odir =~ s/ //g;

###############
###MAIN
###############

if (!-d $odir){
make_path($odir) or die "Error creating directory: $odir";
}

chdir($dir) or die "Could not find directory $dir, $!\n";

unless (-f $file) {
	warn "Could not find input file $file given via -f option";
	pod2usage(-verbose => 0);
}

unless ($type || $file =~ /.fa$/ || $file =~ /.fastq/ || $file =~ /.fasta/) {
	warn "Could not determine file type, please define using the -t option (fasta or fastq)";
	pod2usage(-verbose => 0);
}

open(IN,"<:gzip(autopop)","$file") || die ("Could not open $file!\n@!\n");

my (%embeddings, %kmer, @seqs);
my $x=1; ## Line counter for fastq file
while(<IN>){
    if ($file =~ /\.fastq/ || $type eq "fastq"){
#	print STDERR "Processing fastq file!\n";
        if( (++$x)%2 && $x != 5){
	    chomp (my $raw = $_);
	    push @seqs , $raw;
	}
	else{
#	    print STDERR $x,"\n";
	    $x = 1 if ($x >= 5);
	}
#	$embeddings{"$line[1] $line[2]"}= $line[5];
    }
    elsif ($file =~ /\.fa$/ || $file =~ /\.fasta$/ || $type eq "fasta"){
#	print STDERR "Processing fasta file!\n";
	next if ($_=~m/^>/);
	chomp (my $raw = $_);
	push @seqs, $raw;
#	$embeddings{"seq"}= $raw;
    }
    else{
	warn "Could not determine file type, please define using the -t option (fasta or fastq)";
	pod2usage(-verbose => 0);
    }
}
close(IN);

###Investigate kmer_enrichment
%kmer = %{kmer_enrichment(\@seqs, $klength)};
#foreach my $key (keys %embeddings){
#  %kmer = %{kmer_enrichment($embeddings{$key}, $klength)};
#  %kmer = %{kmer_enrichment(\@seqs, $klength)};
#}

my $total = sum values %kmer;
chdir($odir) or die "Could not find directory $odir, $!\n";

### Print Output
open(KMER,">","$klength\_mers") or die "Could not open file $!\n";
print KMER "$klength\-mer\tCount\tRatio\n";
print KMER "TOTAL\t$total\t1\n";
foreach my $key  (sort {$kmer{$b} <=> $kmer{$a} } keys %kmer) {
    my $ratio = $kmer{$key}/$total;
    print KMER "$key\t$kmer{$key}\t$ratio\n";
}
close(KMER);

chdir($dir) or die "Could not find directory $dir, $!\n";
###############
###POD
###############

__END__

=head1 NAME

kmer_analysis.pl - Simple k-mer count analysis of fasta or fastq files

=head1 SYNOPSIS

kmer_analysis.pl [-f I<FILE>] [-d I<STRING>] [-o I<STRING>] [-k
I<INTEGER>] [-t I<STRING>] [options]

=head1 DESCRIPTION

This program counts k-mers of user defined length in fasta or fastq files.

=head1 OPTIONS

=over 

=item B<-f>

File for processing

=item B<-d>

Working directory

=item B<-o>

Output directory

=item B<-k>

Kmer length to search

=item B<-t>

File type, can either be fasta or fastq

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back


=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################

