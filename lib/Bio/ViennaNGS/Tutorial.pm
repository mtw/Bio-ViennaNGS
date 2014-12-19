# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:41:39 mtw>

package Bio::ViennaNGS::Tutorial;

use 5.12.0;
use Exporter;
use version; our $VERSION = qv('0.12_07');
use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT = ();

our @EXPORT_OK = qw ();

1;
__END__

=head1 NAME

    Tutorial_pipeline01.pl - An example pipeline for the ViennaNGS toolbox

=head2 SYNOPSIS
 
   ./Pipeline.pl path/to/R/libraries
    Where path/to/R/libraries should point to the directory containing ggplot2

=head2 Pipeline 

    This script is a showcase for the usage of ViennaNGS in a real NGS example.
    We start from a file containing ENSEMBL annotation information for human protein-coding genes.
    We are insterested in finding sequence motifs in close proximity to the gene start (50nt upstream, 10nt into the gene).

    The first step is to initialize some variables and generate a chromosome_sizes hash.
    
    C<my $bed	 = 'hg19_highlyexpressed.bed';>
    C<my $name     = (split(/\./,$bed))[0];>
    C<my $upstream = 50;>
    C<my $into     = 10;>
    C<my $outfile  = "$name.ext$upstream\_fromStart_$into\_downstream.bed";>
    C<my $outfile2 = "$name.ext$upstream\_upstream.bed";>
    C<my %sizes = %{fetch_chrom_sizes('hg19')};
    
=cut

=head2 Generate a Bio::ViennaNGS::FeatureChain object

    The bed file of interest is parsed, a feature array is generated and passed on to Bio::ViennaNGS::FeatureChain, which creates a new Moose Object of type FeatureChain, containing the original bed entries
    C<my @featurelist = @{parse_bed6($bed)};
    ### Now we create a Bio::ViennaNGS::FeatureChain from the featurelist above
    my $chain = Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);>
    
=cut

=head2 Extend the existing chain for motif analysis

    The newly created FeatureChain object will now be extended 50nt upstream of the gene start and 10nt into the gene, to retrieve a bed file which contains the putative sequence motifs.
    
    C<my $extended_chain = extend_chain(\%sizes,$chain,0,$into,$upstream,0);>
    
    For later purposes we also extend the whole U6 gene span 50nt upstream.
    C<my $extended_chain2 = extend_chain(\%sizes,$chain,$upstream,0,0,0);>

=cut

=head2 Print extended Bio::ViennaNGS::FeatureChain objects to files

    Extended chains are now print out to make them available for external tools like bedtools.
    C<my $out = $extended_chain->print_chain();
    print $Out $out;
    $out = $extended_chain2->print_chain();
    print $Out2 $out;>

=cut

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

=head2 Sequence analysis

    We now generate FASTA files from the extended bed files using bedtools getfasta method.
    C<my $bedtools = `bedtools getfasta -s -fi hg19_chromchecked.fa -bed $outfile -fo $name.ext$upstream\_fromStart_$into\_downstream.fa`;>
    C<print STDERR "$bedtools\n" if $?;>
    C<$bedtools = `bedtools getfasta -s -fi hg19_chromchecked.fa -bed $outfile2 -fo $name.ext$upstream\_upstream.fa`;>
    C<print STDERR "$bedtools\n" if $?;>

    To analyze putative sequence motifs in the newly generated Fasta files, we use two approaches.
    First we analyze the k-mer content with the Bio::ViennaNGS(kmer_enrichment) method for k-mers of length 6 to 8 nt.

    C<open(IN,"<","$name.ext$upstream\_fromStart_$into\_downstream.fa") || die ("Could not open $name.ext$upstream\_fromStart_$into\_downstream.fa!\n@!\n");>

    C<my @fastaseqs;>
    C<while(<IN>){
          chomp (my $raw = $_);
          next if ($_ =~ /^>/);
          push @fastaseqs, $raw;
    }>
    C<close(IN);>

    C<for (6..8){
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
    }>
=cut

=head1 

    In a second approach we run MEME to retrieve the 20 most over-represented motifs of length 8.
    C<meme hg19_highexpressed.ext50_fromStart_10_downstream.fa -oc MEME_hg19_highexpressed.ext50_fromStart_10_downstream.fa -w 8 -dna -maxsize 1000000000 -nmotifs 20>

    Once the meme run is done, we want to have a nice figure which shows the e-value and site coverage of the top 10 motifs

    C<`./scripts/MEME_xml_motif_extractor.pl -f Example_Pipeline_meme.xml -r $RLIBPATH -t Example_Pipeline`>

=cut

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut

##################################END################################
