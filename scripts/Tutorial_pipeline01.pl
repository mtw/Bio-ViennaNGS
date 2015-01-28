#!/usr/bin/env perl
# Last changed Time-stamp: <2015-01-28 10:38:48 fall>
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
=head1 NAME

Tutorial_pipeline01.pl - An example pipeline for the ViennaNGS toolbox

=head2 SYNOPSIS

./Tutorial_Pipeline01.pl path/to/R/libraries
Where path/to/R/libraries should point to the directory containing ggplot2

=head2 DESCRIPTION 

This script is a showcase for the usage of ViennaNGS in a real NGS
example.  We start from a file containing ENSEMBL annotation
information for human protein-coding genes.  We are insterested in
finding sequence motifs in close proximity to the gene start (50nt
upstream, 10nt into the gene).

=head2 PIPELINE

The first step is to initialize some variables and generate a chromosome_sizes hash.

C<<my $bed	 = 'hg19_highlyexpressed.bed';>>
C<<my $name     = (split(/\./,$bed))[0];>>
C<<my $upstream = 50;>>
C<<my $into     = 10;>>
C<<my $outfile  = "$name.ext$upstream\_fromStart_$into\_downstream.bed";>>
C<<my $outfile2 = "$name.ext$upstream\_upstream.bed";>>
C<<my %sizes = %{fetch_chrom_sizes('hg19')};>>

=cut

my $bed	     = 'hg19_highlyexpressed.bed';
my $name     = (split(/\./,$bed))[0];
my $upstream = 50;
my $into     = 10;
my $outfile  = "$name.ext$upstream\_fromStart_$into\_downstream.bed";

my %sizes = %{fetch_chrom_sizes('hg19')}; ### Requires installation of UCSCs fetchChromSizes script or mysql

=head3 Generate a Bio::ViennaNGS::FeatureChain object

The bed file of interest is parsed, a feature array is generated and
passed on to Bio::ViennaNGS::FeatureChain, which creates a new Moose
Object of type FeatureChain, containing the original bed entries

C<<my @featurelist = @{parse_bed6($bed)};
### Now we create a Bio::ViennaNGS::FeatureChain from the featurelist above
my $chain = Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);>>

=cut

my @featurelist = @{parse_bed6($bed)};
my $chain	= Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);

=head3 Extend the existing chain for motif analysis

The newly created FeatureChain object will now be extended 50nt upstream of the gene start and 10nt into the gene, to retrieve a bed file which contains the putative sequence motifs.

C<<my $extended_chain = extend_chain(\%sizes,$chain,0,$into,$upstream,0);>>

=cut

my $extended_chain  = extend_chain(\%sizes,$chain,0,$into,$upstream,0);

=head3 Print extended Bio::ViennaNGS::FeatureChain objects to files

Extended chains are now print out to make them available for external tools like bedtools.

C<<my $out = $extended_chain->print_chain();
print $Out $out;
close ($Out)>>

=cut

open (my $Out, ">",$outfile) or die "$!";

my $out	= $extended_chain->print_chain();
print $Out $out;

close($Out);

=head3 Summary of so far used methods

=over 4

=item fetch_chrom_sizes<as_string>

Returns a chromosome-sizes hash reference for the specified species, e.g. hg19, mm9, mm10, etc.

=item parse_bed6<as_string>

Reads a bed6 file and returns a feature array.

=item Bio::ViennaNGS::FeatureChain->new()<as_method>

Generates a new Bio::ViennaNGS::FeatureChain object from a feature array

=item Bio::ViennaNGS(extend_chain)<as_method>

Extends a Bio::ViennaNGS::FeatureChain object by given constraints

=back

=cut

=head3 Sequence analysis

 We now generate FASTA files from the extended bed files using bedtools getfasta method.
C<<my $bedtools = `bedtools getfasta -s -fi hg19_chromchecked.fa -bed $outfile -fo $name.ext$upstream\_fromStart_$into\_downstream.fa`;>>
C<<print STDERR "$bedtools\n" if $?;>>
C<<$bedtools = `bedtools getfasta -s -fi hg19_chromchecked.fa -bed $outfile2 -fo $name.ext$upstream\_upstream.fa`;>>
C<<print STDERR "$bedtools\n" if $?;>>

To analyze putative sequence motifs in the newly generated Fasta
files, we use two approaches.  First we analyze the k-mer content with
the Bio::ViennaNGS(kmer_enrichment) method for k-mers of length 6 to 8
nt.

C<<open(IN,"<","$name.ext$upstream\_fromStart_$into\_downstream.fa") || die ("Could not open $name.ext$upstream\_fromStart_$into\_downstream.fa!\n@!\n");>>

C<<my @fastaseqs;>>
C<<while(<IN>){
chomp (my $raw = $_);
next if ($_ =~ /^>/);
push @fastaseqs, $raw;
}>>
C<<close(IN);>>

C<<for (6..8){
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
}>>

=cut

my $bedtools = `bedtools getfasta -s -fi hg19_chromchecked.fa -bed $outfile -fo $name.ext$upstream\_fromStart_$into\_downstream.fa`;
print STDERR "$bedtools\n" if $?;

open(IN,"<","$name.ext$upstream\_fromStart_$into\_downstream.fa") || die ("Could not open $name.ext$upstream\_fromStart_$into\_downstream.fa!\n@!\n");

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

=head3 MEME

In a second approach we run MEME to retrieve the 20 most
over-represented motifs of length 8.  C<meme
hg19_highexpressed.ext50_fromStart_10_downstream.fa -oc
MEME_hg19_highexpressed.ext50_fromStart_10_downstream.fa -w 8 -dna
-nmotifs 20>

Once the meme run is done, we want to have a nice figure which shows
the e-value and site coverage of the top 10 motifs

C<<my $cmd = "perl ../scripts/MEME_xml_motif_extractor.pl -f Example_Pipeline_meme.xml -r $RLIBPATH -t Example_Pipeline";
my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(command => $cmd, verbose => 0);

if(!$success){    
my $this_function = (caller(0))[3];
print STDERR "ERROR: MEME_xml_motif_extractor.pl run unsuccessful\n";
print join "", @$full_buf;
unless ($r) {
warn "If you do not provide a valid R-lib-path and ggplot is not found in the standard R path, this pipeline will not be able to parse the MEME xml output.\n";
pod2usage(-verbose => 0);
}
}>>
=cut

my $cmd = "perl ../scripts/MEME_xml_motif_extractor.pl -f Example_Pipeline_meme.xml -r $RLIBPATH -t Example_Pipeline";
my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(command => $cmd, verbose => 0);

if(!$success){    
    my $this_function = (caller(0))[3];
    print STDERR "ERROR: MEME_xml_motif_extractor.pl run unsuccessful\n";
    print join "", @$full_buf;
    unless ($r) {
	warn "If you do not provide a valid R-lib-path and ggplot is not found in the standard R path, this pipeline will not be able to parse the MEME xml output.\n";
	pod2usage(-verbose => 0);
    }
}


=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut

##################################END################################

