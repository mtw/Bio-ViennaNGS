#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-02-09 16:57:28 fabian>

=head1 NAME

Tutorial_pipeline00.pl - Inferring quantitative and qualitative
parameters from an input BAM file.

=head1 SYNOPSIS

  Tutorial_pipeline00.pl 

=head1 DESCRIPTION

This tutorial illustrates how the libraries L<Bio::ViennaNGS::BamStat> and
L<Bio::ViennaNGS::BamStatSummary> can be used to count mapped reads, read
alignments, single-end and paired-end reads and to check and compare
quality features stored in the BAM file. The latter is exemplify by
applied it to deduce the distribution of edit distance for all read
alignments.

Thereby, I would like to point out that this tutorial does not cover
all feature of L<Bio::ViennaNGS::BamStat> and
L<Bio::ViennaNGS::BamStatSummary>. It is merely meant to illustrate the
principles. For more details on L<Bio::ViennaNGS::BamStat> and
L<Bio::ViennaNGS::BamStatSummary> please refer to their documentation,

=head1 INTRODUCTION

In our toy example we aim to examine the count in an exact way
different types of reads and visualize the distribution of the edit
distance between the aligned reads and the reference genome. To this
end we use the same mapped RNA-seq data as examined in the subsequent
Tutorials, which together are meant as an exemplary analysis
pipeline. As usual quality control of the input data has its natural
place in the beginning.

The input data BAM file can be retrieved
L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS/C1R1.bam>). Please
download the file using your browser or a simple bash command C<wget
http://nibiru.tbi.univie.ac.at/ViennaNGS/C1R1.bam >. From here on, the
input file C1R1.bam is assumed to be accessible in your working
directory.

=head1 PROCEDURE

=head2 Include libraries

 
 use Bio::ViennaNGS::BamStat;
 use Bio::ViennaNGS::BamStatSummary;
 use Data::Dumper;


For this tutorial three special libraries are included. First,
C<<Bio::ViennaNGS::BamStat>> provides methods to read key aspects
concerning quality and quantity information from a defined BAM file
and stores its essential information in a data object. Second, to get
an impression of the data stored the standard perl library
C<<Data::Dumper>> is often usefully. Third,
C<<Bio::ViennaNGS::BamStatSummary>> provides methods to compare and
visualize stored in the BamStat object.

=cut

 use strict;
 use warnings;
 use Getopt::Long qw( :config posix_default bundling no_ignore_case );
 use Pod::Usage;
 use Data::Dumper;
 use File::Basename;
 use Path::Class;
 use Bio::ViennaNGS::BamStat;
 use Bio::ViennaNGS::BamStatSummary;


=head2 Define control variables

 @bams     = qw# C1R1.bam #;
 $odir     = '.';

 $edit_control     = 1;
 $segemehl_control = 1;


In the first step of this tutorial, control and parameter variables
are set. The array C<@bams> holds a list to all BAM files intended to
be used in this analysis. We restrict ourself here to one file named
C1R1.bam, which should be accessible in the current working directory.

Since Tutorial_pipeline00.pl produces several output file, with fixed
file names, an output directory has to be specified. This is done in
the C<$odir> variable, setting it here to the current working
directory.  Please not that if files with same names do already exist
in this particular directory, they will be overwritten.

Some methods in L<Bio::ViennaNGS::BamStatSummary> use the
L<Bio::ViennaNGS::BamStatSummary> use the C<Statistics::R>
library. Therefore, a absolute and valid path to the a working version
of R has to be specified. Here, as in the most standard Linux
installations, the path is set to /usr/bin/R.

The next control variable C<$edit_control> flags if in the course of
populating the data object by L<Bio::ViennaNGS::BamStat> information on
the edit distance of the read alignments should be stored. If so it
will be usable to visualize this information in a subsequent step by
L<Bio::ViennaNGS::BamStatSummary>.

L<Bio::ViennaNGS::BamStat> and L<Bio::ViennaNGS::BamStatSummary> are
in principle compatible with any BAM file from any read
aligner. Nevertheless one has to be aware that some aligner differ in
respect to the BAM dialect or with respect to the information stored
in the BAM file. Therefore, we introduce here a special flag to use
the auxiliary information stored in the BAM file produced by
I<segemehl>. Toggle it to '1' if your input file is from this origin,
as it is the case of the provided C1R1.bam. Otherwise set to '0'.

=cut

my $VERBOSE = 0;
my @bams     = qw# C1R1.bam #;
next unless ($bams[0] eq $bams[-1] && $bams[0] eq "C1R1.bam");
my $odir     = './';
my %data     = ();

###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
					   "help|h"    => sub{pod2usage(-verbose => 1)},
					   "man|m"     => sub{pod2usage(-verbose => 2)},
					   "verbose"   => sub{ $VERBOSE++ }
					  );


=head2 Creating new BamStatSummary object.

 $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						   outpath        => $odir,
						   is_segemehl    => $segemehl_control,
						   control_edit   => $edit_control,
						  );


Initialize new BamStatSummary object capable of representing data from
 all I<segemehl> flavored BAM files (C<<< is_segemehl => 1 >>>)
 specified in @bams, setting the output directory to $odir ('./'),
 where beside standard read quantification also the edit distance
 (C<<< control_edit => 1 >>>) of each read will be stored.

=cut

my $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						     outpath        => $odir,
                                                     is_segemehl    => 1,
                                                     control_edit   => 1,
                                                     control_flag   => 1,
                                                     control_score  => 1,

						    );


=head2 Read-in BAM files @bams

In the next step the initialized data object has to be populated with
data. Therefore, each BAM file has to be read an the essential data
has to be extracted. 

 $bamsummary->populate_data();

Thereby, for each BAM file in @bams the method C<< new >> from
L<BIO::ViennaNGS::BamStat> is called like this, 

 $bo = Bio::ViennaNGS::BamStat->new(bam => $bamfile);

This calls internally L<Bio::DB::Sam> library within
L<Bio::ViennaNGS::BamStat>.

To examine the content of this new object, please use C< print
Dumper($bamsummary); >. As you will see its content depends on the way
it is initialized, and the data specified to be stored.

=cut

$bamsummary->populate_data();

=head2 Quantify data from $bamsummary

The most basic step in the course of the assessment of input BAM files
is the quantification. Namely, how many reads are uniquely or multiply
mapped? how many alignments are there? Are there single-end or
paired-end reads, and if the latter how many pairs are complete? To
this end we can compile needed quantitative information for all reads
stored in $bamsummary by,

 $bamsummary->populate_countStat();

Again, you can use C<print Dumper($bamsummary);> to examine the
object.

=cut

$bamsummary->populate_countStat();

=head2 Produce output for read quantification.

In the next step we will out put the compile information into a file
in $odir.

 $bamsummary->dump_countStat("csv");

The output file format is *.csv which can easily be screened with any
text editor or spreadsheet program.

=cut

$bamsummary->dump_countStat("csv");

=head2 Plot read quantification.

Beside the summarized read quantification in tabular form. It can be
useful to plot the numbers. This can help to get a quick overview of
the consistency of different examined samples.

 $bamsummary->make_BarPlot();

This creates a barplot for the read quantification. The file format
is *.pdf, and again the file is created in $odir.

=cut

$bamsummary->make_BarPlot();

=head2 Plot edit distance distribution.

Finally, we would like to gain a quick overview of the quality of
different mapped RNA-seq samples. Therefore we like to plot the
distribution of edit distances for all reads aligned to the reference
genome for all samples in @bams.

 $bamsummary->make_BoxPlot("data_edit") if( $bamsummary->has_control_edit );

=cut


$bamsummary->make_BoxPlot("data_edit")  if( $bamsummary->control_edit  && $bamsummary->has_control_edit );

=head2 Summary

In the previous seven sections we used L<Bio::ViennaNGS::BamStat> and
L<Bio:ViennaNGS::BamStatSummary> to extract, store, summarize, and
visualize quantity and quality data stored in a BAM file. Only
exemplary features of the library were illustrated. It's modular
architecture allows easily to extend its functionality. Further useful
functions are all ready implemented in the corresponding utility
bam_quality_stat.pl. Further can be implemented in a customized manner
according to own needs.

=head1 AUTHOR

Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=cut
