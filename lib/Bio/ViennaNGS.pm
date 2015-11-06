# -*-CPerl-*-
# Last changed Time-stamp: <2015-11-06 13:36:49 mtw>

package Bio::ViennaNGS;

use version; our $VERSION = qv('0.17_01');

1;

=head1 NAME

Bio::ViennaNGS - A Perl distribution for Next-Generation Sequencing
(NGS) data analysis

=head1 DESCRIPTION

Bio::ViennaNGS is a distribution of Perl modules and utilities for
building efficient Next-Generation Sequencing (NGS) analysis
pipelines. It covers various aspects of NGS data analysis, including
(but not limited to) conversion of sequence annotation, evaluation of
mapped data, expression quantification and visualization.

The main Bio::ViennaNGS module is shipped with a complementary set of
(sub)modules:

=over

=item L<Bio::ViennaNGS::AnnoC>: A Moose interface for storage and
conversion of sequence annotation data.

=item L<Bio::ViennaNGS::Bam>: Routines for high-level manipulation of
BAM files.

=item L<Bio::ViennaNGS::BamStat>: A L<Moose> based class for
collecting mapping statistics.

=item L<Bio::ViennaNGS::BamStatSummary>: A L<Moose> interface for
processing L<Bio::ViennaNGS::BamStat> objects on a set of BAM files.

=item L<Bio::ViennaNGS::Bed>: A L<Moose> interface for manipulation of
genomic interval data in BED format.

=item L<Bio::ViennaNGS::BedGraphEntry>: A L<Moose> interface for
storing genomic cvoerage and interval data in bedGraph format.

=item L<Bio::ViennaNGS::Expression>: An object oriented interface for
read-count based gene expression analysis.

=item L<Bio::ViennaNGS::ExtFeature>: A L<Moose> wrapper for extended
BED6 entries.

=item L<Bio::ViennaNGS::Fasta>: Routines for accessing genomic
sequences implemented through a L<Moose> interface to
L<Bio::DB::Fasta>.

=item L<Bio::ViennaNGS::Feature>: A L<Moose> based BED6 wrapper.

=item L<Bio::ViennaNGS::FeatureChain>: Yet another L<Moose> class for
chaining gene annotation features.

=item L<Bio::ViennaNGS::FeatureInterval>: A L<Moose> interface for
handling elementary genomic intervals, corresponding to BED3.

=item L<Bio::ViennaNGS::FeatureIO>: A L<Moose> interface for efficient
input/output handling of genomic annotation formats.

=item L<Bio::ViennaNGS::FeatureLine>: An abstract L<Moose> class for
combining several L<Bio::ViennaNGS::FeatureChain> objects.

=item L<Bio::ViennaNGS::MinimalFeature>: A L<Moose> interface for
handling elementary gene annotation, corresponding to BED4.

=item L<Bio::ViennaNGS::Peak>: A L<Moose> interface for identification
and characterization of peaks/enriched regions in RNA-seq data.

=item L<Bio::ViennaNGS::SpliceJunc>: A collection of routines for
alternative splicing analysis.

=item L<Bio::ViennaNGS::Tutorial>: A comprehensive tutorial of the
  L<Bio::ViennaNGS> core routines with real-world NGS data.

=item L<Bio::ViennaNGS::UCSC>: Routines for visualization of genomics
data with the UCSC genome browser.

=item L<Bio::ViennaNGS::Util>: A collection of wrapper routines for
commonly used third-party NGS utilities, code for normalization of
gene expression values based on read count data and a set of utility
functions.

=back

=head1 UTILITIES

L<Bio::ViennaNGS> comes with a collection of command line utilities
for accomplishing routine tasks often required in NGS data
processing. These utilities serve as reference implementation of the
routines implemented throughout the modules and can readily be used
for atomic tasks in NGS data processing:

=over

=item F<assembly_hub_constructor.pl>: The UCSC genome browser offers
the possibility to visualize any organism (including organisms that
are not included in the standard UCSC browser bundle) through hso
called 'Assembly Hubs'. This script constructs Assembly Hubs from
genomic sequence and annotation data.

=item F<bam_split.pl>: Split (paired-end and single-end) BAM alignment
files by strand and compute statistics. Optionally create BED output,
as well as normalized bedGraph and bigWig files for coverage
visualization in genome browsers (see dependencies on third-patry
tools below).

=item F<bam_to_bigWig.pl>: Produce bigWig coverage profiles from
(aligned) BAM files, explicitly considering strandedness. The most
natural use case of this tool is to create strand-aware coverage
profiles in bigWig format for genome browser visualization.

=item F<bam_uniq.pl>: Extract unique and multi mapping reads from BAM
alignment files and create a separate BAM file for both uniqe (.uniq.)
and multi (.mult.) mappers.

=item F<bed2bedGraph.pl>: Convert BED files to (strand specific)
bedGraph files, allowing additional annotation and automatic
generation of bedGraph files which can easily be converted to big-type
files for easy UCSC visualization.

=item F<bed2nt2aa.pl>: Provide nucleotide and amino acid sequences for
BED6 intervals.

=item F<extend_bed.pl>: Extend genomic features in BED files by a
certain number of nucleotides, either on both sides or specifically at
the 5' or 3' end, respectively.

=item F<gff2bed.pl>: Convert RefSeq GFF3 annotation files to BED12
format. Individual BED12 files are created for each feature type
(CDS/tRNA/rRNA/etc.). Tested with RefSeq bacterial GFF3 annotation.

=item F<kmer_analysis.pl>: Count k-mers of predefined length in FastQ
and Fasta files

=item F<MEME_XML_motif_extractor.pl>: Compute simple statistics from
MEME XML output and return a list of found motifs with the number of
sequences containing those motifs as well as nice ggplot graphs.

=item F<newUCSCdb.pl>: Create a new genome database to a locally
installed instance of the UCSC genome browser in order to add a novel
organism for visualization. Based on L<this Genomewiki
article|http://genomewiki.ucsc.edu/index.php/Building_a_new_genome_database>.

=item F<normalize_multicov.pl>: Compute normalized expression data in
RPKM and TPM from (raw) read counts in bedtools multicov format. TPM
reference: Wagner et al, Theory Biosci. 131(4), pp 281-85 (2012).

=item F<rnaseq_peakfinder.pl>: Find and characterize peaks/enriched
regions of certain size and coverage in RNA-seq data.

=item F<sj_visualizer.pl>: Convert splice junctions from mapped
RNA-seq data in segemehl BED6 splice junction format to BED12 for easy
visualization in genome Browsers.

=item F<splice_site_summary.pl>: Identify and characterize splice
junctions from RNA-seq data by intersecting them with annotated splice
junctions.

=item F<trim_fastq.pl>: Trim sequence and quality string fields in a
Fastq file by user defined length.

=back

=head1 DEPENDENCIES

The L<Bio::ViennaNGS> modules and classes depend on a set of Perl
modules, some of which are part of the Perl core distribution:

=over

=item L<Bio::Perl> >= 1.00690001

=item L<Bio::DB::Sam> >= 1.37

=item L<Bio::DB::Fasta>

=item L<Bio::Tools::GFF>

=item L<File::Basename>

=item L<File::Share>

=item L<File::Slurp>

=item L<File::Temp>

=item L<List::Util>

=item L<Path::Class>

=item L<IPC::Cmd>

=item L<Carp>

=item L<Template>

=item L<Moose>

=item L<Moose::Util::TypeConstraints>

=item L<namespace::autoclean>

=item L<MooseX::Clone>

=item L<MooseX::InstanceTracking>

=item L<Tie::Hash::Indexed>

=item L<Test::Files>

=item L<Test::File::Contents>

=back

In addition the following modules are required by the L<Bio::ViennaNGS> utilities:

=over

=item L<PerlIO::gzip>

=item L<Math::Round>

=item L<XML::Simple>

=item L<Statistics::R>

=back

L<Bio::ViennaNGS> depends on a set of third-party tools and libraries which
are required for specific filtering and file format conversion tasks as
well as for building internally used Perl modules:

=over

=item L<bedtools2|https://github.com/arq5x/bedtools2>

=item F<bedGraphToBigWig>, F<fetchChromSizes>, F<faToTwoBit> from the
  L<UCSC Genome Browser
  applications|http://hgdownload.cse.ucsc.edu/admin/exe>

=item L<R|http://www.r-project.org/>

=item samtools E<lt>=v0.1.19 for building L<Bio::DB::Sam>.

=back

Please ensure that all third-party utilities are available on your
system, and that hey can be found and executed by the Perl
interpreter.

=head1 SOURCE AVAILABILITY

Source code for this distribution is available from the L<ViennaNGS
Github repository|https://github.com/mtw/Bio-ViennaNGS>.

=head1 PAPERS

If the L<Bio::ViennaNGS> suite is useful for your work or if you use
any component of the distribution in a custom pipeline, please cite
the following publication:

B<"ViennaNGS - A toolbox for building efficient next-generation sequencing
analysis pipelines">

I<Michael T. Wolfinger, Joerg Fallmann, Florian Eggenhofer and Fabian Amman>

F1000Research 2015, 4:50 (doi: L<10.12688E<sol>f1000research.6157.2|http://dx.doi.org/10.12688/f1000research.6157.2>)

=head1 NOTES

The L<Bio::ViennaNGS> suite is actively developed and tested on
different flavours of Linux and Mac OS X. We have taken care of
writing platform-independent code that should run out of the box on
most UNIX-based systems, however we do not have access to machines
running Microsoft Windows. As such, we have not tested and will not
test Windows compatibility.

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS::AnnoC>

=item L<Bio::ViennaNGS::Bam>

=item L<Bio::ViennaNGS::BamStat>

=item L<Bio::ViennaNGS::BamStatSummary>

=item L<Bio::ViennaNGS::Bed>

=item L<Bio::ViennaNGS::BedGraphEntry>

=item L<Bio::ViennaNGS::Expression>

=item L<Bio::ViennaNGS::ExtFeature>

=item L<Bio::ViennaNGS::Fasta>

=item L<Bio::ViennaNGS::Feature>

=item L<Bio::ViennaNGS::FeatureChain>

=item L<Bio::ViennaNGS::FeatureIO>

=item L<Bio::ViennaNGS::FeatureInterval>

=item L<Bio::ViennaNGS::FeatureLine>

=item L<Bio::ViennaNGS::MinimalFeature>

=item L<Bio::ViennaNGS::Peak>

=item L<Bio::ViennaNGS::SpliceJunc>

=item L<Bio::ViennaNGS::Tutorial>

=item L<Bio::ViennaNGS::UCSC>

=item L<Bio::ViennaNGS::Util>

=back

=head1 AUTHORS

=over

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=item Joerg Fallmann E<lt>fall@tbi.univie.ac.atE<gt>

=item Florian Eggenhofer E<lt>florian.eggenhofer@tbi.univie.ac.atE<gt>

=item Fabian Amman E<lt>fabian@tbi.univie.ac.at<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2015 Michael T. Wolfinger
E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut


