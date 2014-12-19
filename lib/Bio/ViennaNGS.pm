# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:30:22 mtw>

package Bio::ViennaNGS;

use version; our $VERSION = qv('0.12_07');

1;

=head1 NAME

Bio::ViennaNGS - A Perl distribution for Next-Generation Sequencing
(NGS) data analysis

=head1 DESCRIPTION

Bio::ViennaNGS is a distribution of Perl modules and utilities for building
efficient Next-Generation Sequencing (NGS) analysis pipelines. It covers
various aspects of NGS data analysis, including (but not limited to)
conversion of sequence annotation, evaluation of mapped data, expression
quantification and visualization.

The main Bio::ViennaNGS module is shipped with a complementary set of 
(sub)modules:

=over 

=item L<Bio::ViennaNGS::Fasta>: Routines for accessing genomic
sequences implemented through a Moose interface to L<Bio::DB::Fasta>.

=item L<Bio::ViennaNGS::AnnoC>: A Moose interface for storage and
  conversion of sequence annotation data.

=item L<Bio::ViennaNGS::SpliceJunc>: A collection of routines for
  alternative splicing analysis.

=item L<Bio::ViennaNGS::UCSC>: Routines for visualization of genomics
data with the UCSC genome browser.

=item L<Bio::ViennaNGS::Util>: Commodity routines used throughout the
    L<Bio::ViennaNGS> modules and utilities

=item L<Bio::ViennaNGS::MinimalFeature>: A Moose interface for
handling elementary gene annotation.

=item L<Bio::ViennaNGS::Feature>: A Moose-ish BED6 wrapper.  

=item L<Bio::ViennaNGS::FeatureChain>: Yet another Moose class for
chaining gene annotation features.

=back

L<Bio::ViennaNGS> comes with a collection of command line utilities
for accomplishing routine tasks often required in NGS data
processing. These utilities serve as reference implementation of the
routines implemented throughout the modules and can readily be used
for atomic tasks in NGS data processing:

=over

assembly_hub_constructor.pl:
The UCSC genome browser offers the possibility to visualize any
organism (including organisms that are not included in the standard
UCSC browser bundle) through hso called 'Assembly Hubs'. This script
constructs Assembly Hubs from genomic sequence and annotation data.

bam_split.pl: 
Split (paired-end and single-end) BAM alignment files by strand and compute
statistics. Optionally create BED output, as well as normalized bedGraph
and bigWig files for coverage visualization in genome browsers (see
dependencies on third-patry tools below).

bam_to_bigWig.pl:
Produce bigWig coverage profiles from (aligned) BAM files, explicitly
considering strandedness. The most natural use case of this tool is to
create strand-aware coverage profiles in bigWig format for genome browser
visualization.

bam_uniq.pl:
Extract unique and multi mapping reads from BAM alignment files and create
a separate BAM file for both uniqe (.uniq.) and multi (.mult.) mappers.

bed2bedGraph.pl:
Convert BED files to (strand specific) bedGraph files, allowing
additional annotation and automatic generation of bedGraph files which can
easily be converted to big-type files for easy UCSC visualization.

extend_bed.pl:
Extend genomic features in BED files by a certain number of nucleotides,
either on both sides or specifically at the 5' or 3' end, respectively.

gff2bed.pl:
Convert RefSeq GFF3 annotation files to BED12 format. Individual BED12
files are created for each feature type (CDS/tRNA/rRNA/etc.). Tested with
RefSeq bacterial GFF3 annotation.  

kmer_analysis.pl:
Count k-mers of predefined length in FastQ and Fasta files

MEME_XML_motif_extractor.pl:
Compute simple statistics from MEME XML output and return a list of found
motifs with the number of sequences containing those motifs as well as nice
ggplot graphs.

motiffinda.pl: 
Find motifs in annotated sequence features. The motif can be provided as
regular expression.

newUCSCdb.pl:
Create a new genome database (ie. add a novel organism) in a local
instance of the UCSC genome browser.

normalize_multicov.pl:
Compute normalized expression data in TPM/RPKM from (raw) read counts in
bedtools multicov format. TPM reference: Wagner et al, Theory
Biosci. 131(4), pp 281-85 (2012)

sj_visualizer.pl:
Convert splice junctions from mapped RNA-seq data in segemehl BED6 splice
junction format to BED12 for easy visualization in genome Browsers.

splice_site_summary.pl:
Identify and characterize splice junctions from RNA-seq data by
intersecting them with annotated splice junctions.

newUCSCdb.pl:
Create a new genome database for a locally installed instance of the UCSC
genome browser. Based on
http://genomewiki.ucsc.edu/index.php/Building_a_new_genome_database

trim_fastq.pl:
Trim sequence and quality string fields in a Fastq file
by user defined length.

=back

=head1 DEPENDENCIES

=over 7

=item  L<Bio::Perl> >= 1.00690001

=item  L<BIO::DB::Sam> >= 1.39

=item  L<File::Basename>

=item  L<File::Temp>

=item  L<Path::Class>

=item  L<IPC::Cmd>

=item  L<Carp>

=back

L<Bio::ViennaNGS> uses third-party tools for computing intersections
of BED files: F<bedtools intersect> from the
L<BEDtools|http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html>
suite is used to compute overlaps and F<bedtools sort> is used to sort
BED output files. Make sure that those third-party utilities are
available on your system, and that hey can be found and executed by
the Perl interpreter. We recommend installing the latest version of
L<BEDtools|https://github.com/arq5x/bedtools2> on your system.

=head1 NOTES

The L<Bio::ViennaNGS> suite is actively developed and tested on
different flavours of Linux and Mac OS X. We have taken care of
writing platform-independent code that should run out of the box on
most UNIX-based systems, however we do not have access to machines
running Microsoft Windows. To put it straight: We have not tested and
will not test Windows compatibility. 

=head1 SEE ALSO

=over 8

=item L<Bio::ViennaNGS::Util>

=item L<Bio::ViennaNGS::AnnoC>

=item L<Bio::ViennaNGS::UCSC>

=item L<Bio::ViennaNGS::SpliceJunc>

=item L<Bio::ViennaNGS::Fasta>

=item L<Bio::ViennaNGS::MinimalFeature>

=item L<Bio::ViennaNGS::Feature>

=item L<Bio::ViennaNGS::FeatureChain>

=back

=head1 AUTHORS

=over

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=item Jörg Fallmann E<lt>fall@tbi.univie.ac.atE<gt>

=item Florian Eggenhofer E<lt>florian.eggenhofer@tbi.univie.ac.atE<gt>

=item Fabian Amman E<lt>fabian@tbi.univie.ac.at<gt>


=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut


