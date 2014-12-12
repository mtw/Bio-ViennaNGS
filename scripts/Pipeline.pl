#!/usr/bin/env perl
# Last changed Time-stamp: <2014-12-12 16:45:22 fall>
# AUTHOR: Joerg Fallmann <joerg.fallmann@univie.ac.at>

###############
###Use stuff
###############
use strict;
use warnings;
use Cwd;
use Pod::Usage;
use File::Path qw(make_path remove_tree);
use File::Basename;
use Math::Round;
use Bio::ViennaNGS qw(parse_bed6 extend_chain kmer_enrichment);
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::FeatureChain;
use List::Util qw(sum);
use Data::Dumper;

###############
### MAIN
###############

my $RLIBPATH = shift;

=head1 NAME

    Pipeline.pl - An example pipeline for the ViennaNGS toolbox

=head1 SYNOPSIS
 
   ./Pipeline.pl path/to/R/libraries

=head1 Pipeline

    This script is a showcase for the usage of ViennaNGS in a real NGS example.
    We start from a file containing ENSEMBL annotation information for U6-RNA coding genes.
    We are insterested in finding sequence motifs in close proximity to the gene start (200nt upstream), which mark annotated genes as transcribed by Polymerase III.

    First step is to initialize some variables we need, and generate a chromosome_sizes hash.
    
    C<my $genome = 'mm9.chrom.sizes';
    my $bed = 'U6.bed';
    my $upstream = 200;
    my $downstream = 20;
    my $outfile = 'U6.ext200up_20in.bed';

    open (my $Genome, "<:gzip(autopop)",$g) or die "$!";
    open (my $Out, ">",$o) or die "$!";

    my %sizes;

    while (<$Genome>){
        chomp $_;
        my ($chr,$size)=split (/\t/,$_);
       $sizes{$chr}=$size;
   }>
    
=cut

my $genome   = 'hg19.chrom.size';
my $bed	     = 'U6.bed';
my $upstream = 200;
my $into     = 20;
my $outfile  = 'U6.ext200fromStart_20downstream.bed';
my $outfile2 = 'U6.ext200upstream.bed';

open (my $Genome, "<:gzip(autopop)",$genome) or die "$!";
open (my $Out, ">",$outfile) or die "$!";
open (my $Out2, ">",$outfile2) or die "$!";

my %sizes;

while (<$Genome>){
    chomp $_;
    my ($chr,$size) = split (/\t/,$_);
    $sizes{$chr}    = $size;
}

#print Dumper(\%sizes);
=head1 Pipeline2

    Then the bed file of interest is parsed, features are generated and passed on to Bio::ViennaNGS::FeatureChain, which creates a new Moose Object of type FeatureChain, containing the original bed entries
    C<my @featurelist = @{parse_bed6($bed)};
    ### Now we create a Bio::ViennaNGS::FeatureChain from the featurelist above
    my $chain = Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);>
    
=cut

my @featurelist = @{parse_bed6($bed)};
my $chain	= Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);

#print Dumper($chain);
=head1 Pipeline3

    The newly created FeatureChain object will now be extended 200nt upstream of the gene start and 20nt into the gene, to retrieve a bed file which contains the putative sequence motifs.
    
    C<my $extended_chain = extend_chain(\%sizes,$chain,0,$into,$upstream,0);>
    
    For later purposes we also extend the whole U6 gene span 200nt upstream.
    C<my $extended_chain2 = extend_chain(\%sizes,$chain,200,0,0,0);>

=cut

my $extended_chain  = extend_chain(\%sizes,$chain,0,$into,$upstream,0);
my $extended_chain2 = extend_chain(\%sizes,$chain,200,0,0,0);

=head1 Pipeline4

    Extended chains are now print out to make them available for external tools like bedtools.
    C<my $out = $extended_chain->print_chain();
    print $Out $out;
    $out = $extended_chain2->print_chain();
    print $Out2 $out;>

=cut

my $out	= $extended_chain->print_chain();
print $Out $out;
$out	= $extended_chain2->print_chain();
print $Out2 $out;

=head2 Methods

=over 4

=item parse_bed6<as_string>

    Reads a bed6 file and returns a feature array.

=item Bio::ViennaNGS::FeatureChain->new()<as_method>

    Generates a new Bio::ViennaNGS::FeatureChain object from a feature array

=item Bio::ViennaNGS(extend_chain)<as_method>

    Extends a Bio::ViennaNGS::FeatureChain object by given constraints

=back

=cut

=head1 Pipeline5

    We now generate FASTA files from the extended bed files using bedtools getfasta method.
    I<`bedtools getfasta -fi hg19_chromchecked.fa -bed U6.ext200fromStart_20downstream.bed -fo U6.ext200fromStart_20downstream.fa -s`
    `bedtools getfasta -fi hg19_chromchecked.fa -bed U6.ext200upstream.bed -fo U6.ext200upstream.fa -s`>

    To analyze the sequence motif content of the newly generated Fasta files, we use two approaches.
    First we analyze the k-mer content with the Bio::ViennaNGS(kmer_enrichment) method for k-mers of length 7 to 10 nt.

    C<open(IN,"<","U6.ext200fromStart_20downstream.fa") || die ("Could not open $file!\n@!\n");

    my @fastaseqs;
    while(<IN>){
        chomp (my $raw = $_);
	push @fastaseqs, $raw;
    }
    close(IN);

    for (7..10){
        %kmer = %{kmer_enrichment(\@seqs, $_)};
        my $total = sum values %kmer;
        ### Print Output
        open(KMER,">","$_\_mers") or die "Could not open file $!\n";
        print KMER "$_\-mer\tCount\tRatio\n";
        print KMER "TOTAL\t$total\t1\n";
        foreach my $key  (sort {$kmer{$b} <=> $kmer{$a} } keys %kmer) {
    	    my $ratio = $kmer{$key}/$total;
    	    print KMER "$key\t$kmer{$key}\t$ratio\n";
        }
        close(KMER);
    }>

=cut
my $bedtools = `bedtools getfasta -fi hg19_chromchecked.fa -bed U6.ext200fromStart_20downstream.bed -fo U6.ext200fromStart_20downstream.fa -s`;
print STDERR "$bedtools\n" if $?;
$bedtools = `bedtools getfasta -fi hg19_chromchecked.fa -bed U6.ext200upstream.bed -fo U6.ext200upstream.fa -s`;
print STDERR "$bedtools\n" if $?;

open(IN,"<","U6.ext200fromStart_20downstream.fa") || die ("Could not open U6.ext200fromStart_20downstream.fa!\n@!\n");

my @fastaseqs;
while(<IN>){
    chomp (my $raw = $_);
    next if ($_ =~ /^>/);
    push @fastaseqs, $raw;
}
close(IN);

for (7..10){
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

=head1 Pipeline6

    In a second approach we run MEME to retrieve the 20 most over-represented motifs of length 7-10.
    meme.bin U6.ext200fromStart_20downstream.fa -oc Meme_U6 -minw 7 -maxw 10 -dna -nmotifs 20 -maxsize 10000000

    Once the meme run is done, we want to have a nice figure which shows the e-value and site coverage of the top 10 motifs

    C<'./scripts/MEME_xml_motif_extractor.pl -f  meme.xml -r $RLIBPATH -t U6'>

=cut

`perl ../scripts/MEME_xml_motif_extractor.pl -f meme.xml -r \$RLIBPATH -t U6`;

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################

