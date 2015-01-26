#!/usr/bin/env perl
# Last changed Time-stamp: <2015-01-26 19:09:13 fall>
# AUTHOR: Joerg Fallmann <joerg.fallmann@univie.ac.at>

###############
###Use stuff
###############
use strict;
use warnings;
use Cwd;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use File::Path qw(make_path remove_tree);
use Math::Round;
use Bio::ViennaNGS::Util qw(parse_bed6 extend_chain kmer_enrichment fetch_chrom_sizes);
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::FeatureChain;
use List::Util qw(sum);
use Data::Dumper;
use IPC::Cmd qw(can_run run);
###############
###Variables
###############

my $VERBOSE = 0;
my ($r,$RLIBPATH)=('','');

###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
    "rpath|r=s" => \$RLIBPATH,
    "help|h"    => sub{pod2usage(-verbose => 1)},
    "man|m"     => sub{pod2usage(-verbose => 2)},      
    "verbose"   => sub{ $VERBOSE++ }
    );

###############
### MAIN
###############

my $bed	     = 'hg19_highlyexpressed.bed';
my $name     = (split(/\./,$bed))[0];
my $upstream = 50;
my $outfile = "$name.ext$upstream\_upstream.bed";

my %sizes = %{fetch_chrom_sizes('hg19')}; ### Requires installation of UCSCs fetchChromSizes script or mysql

my @featurelist = @{parse_bed6($bed)};
my $chain	= Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);

my $extended_chain = extend_chain(\%sizes,$chain,$upstream,0,0,0);

open (my $Out, ">",$outfile) or die "$!";

my $out = $extended_chain->print_chain();
print $Out $out;

close($Out);

print STDERR $extended_chain->type();

my $bedtools = `bedtools getfasta -s -fi hg19_chromchecked.fa -bed $outfile -fo $name.ext$upstream\_upstream.fa`;
print STDERR "$bedtools\n" if $?;

open(IN,"<","$name.ext$upstream\_upstream.fa") || die ("Could not open $name.ext$upstream\_upstream.fa!\n@!\n");

my @fastaseqs;
while(<IN>){
    chomp (my $raw = $_);
    next if ($_ =~ /^>/);
    push @fastaseqs, $raw;
}
close(IN);

for (6..8){
    my %kmer = %{kmer_enrichment(\@fastaseqs, $_)};
    my $total = sum values %kmer;
    ### Print Output
    open(KMER,">","$_\_mers") or die "Could not open file $_\_mers$!\n";
    print KMER "$_\-mer\tCount\tRatio\n";
    print KMER "TOTAL\t$total\t1\n";
    foreach my $key  (sort {$kmer{$b} <=> $kmer{$a} } keys %kmer) {
	my $ratio = nearest(.0001,$kmer{$key}/$total);
	print KMER "$key\t$kmer{$key}\t$ratio\n";
    }
    close(KMER);
}

#my $cmd = "perl ../scripts/MEME_xml_motif_extractor.pl -f Example_Pipeline_meme.xml -r $RLIBPATH -t Example_Pipeline";
#my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(command => $cmd, verbose => 0);
#
#if(!$success){    
#    my $this_function = (caller(0))[3];
#    print STDERR "ERROR: MEME_xml_motif_extractor.pl run unsuccessful\n";
#    print join "", @$full_buf;
#    unless ($r) {
#	warn "If you do not provide a valid R-lib-path and ggplot is not found in the standard R path, this pipeline will not be able to parse the MEME xml output.\n";
#	pod2usage(-verbose => 0);
#    }
#}

###############
###POD
###############

__END__

=head1 NAME

Tutorial_Pipeline02.pl - Another example pipeline for the ViennaNGS toolbox

=head1 SYNOPSIS
    
./Tutorial_Pipeline02.pl

=head1 DESCRIPTION 

This script is a showcase for the usage of ViennaNGS in a real NGS example.
We start from a file containing ENSEMBL annotation information for human protein-coding genes, which have a read pileup of at least 1001 reads in an ENCODE dataset mapped with segemehl.
We are insterested in generating a UCSC bigwig track for those genes and the region 50nt upstream of the gene start.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################

