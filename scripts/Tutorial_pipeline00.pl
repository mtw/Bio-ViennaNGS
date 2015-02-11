#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-02-11 17:36:57 mtw>

=head1 NAME

Tutorial_pipeline00.pl - Inferring quantitative and qualitative
parameters from an input BAM file.

=head1 SYNOPSIS

  Tutorial_pipeline00.pl

=head1 DESCRIPTION

This tutorial illustrates how the libraries L<Bio::ViennaNGS::BamStat>
and L<Bio::ViennaNGS::BamStatSummary> can be used to count mapped
reads, read alignments, single-end and paired-end reads and to check
and compare quality features stored in the BAM file. The latter is
exemplified by using the module to deduce the distribution of edit
distances for all read alignments.

However, I would like to point out that this tutorial does not cover
all features of L<Bio::ViennaNGS::BamStat> and
L<Bio::ViennaNGS::BamStatSummary>. It is merely meant to illustrate
the basic usage principles. For more details please refer to the
documentation of L<Bio::ViennaNGS::BamStat> and
L<Bio::ViennaNGS::BamStatSummary>.

=head1 INTRODUCTION

In this tutorial we examine the count of different types of reads in
an exact way and visualize the distribution of the edit distances
between the aligned reads and the reference genome. 

The input data BAM file can be downloaded
L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS/C1R1.bam>). From here
on, the input file C1R1.bam is assumed to be accessible in your
working directory.

=head1 PROCEDURE

=head2 Include libraries
 
 use Bio::ViennaNGS::BamStat;
 use Bio::ViennaNGS::BamStatSummary;

For this tutorial following C<<Bio::ViennaNGS>> libraries are
included. C<<Bio::ViennaNGS::BamStat>> provides methods to read key
aspects concerning quality and quantity information from a defined BAM
file and stores essential information in a data
object. C<<Bio::ViennaNGS::BamStatSummary>> provides methods to
compare and visualize data stored in the BamStat object.

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
are set. The array C<@bams> holds a list of all BAM files to be
analyzed. This tutorial is restricted to a single file C1R1.bam, which
should be accessible in the current working directory.

Since Tutorial_pipeline00.pl produces several output file, with fixed
file names, an output directory has to be specified. This is done in
the C<$odir> variable, per default set to the current working
directory.  Please not that files with same name in this directory
will be overwritten.

Some methods in L<Bio::ViennaNGS::BamStatSummary> use the
C<Statistics::R> library. Therefore, a valid path to a working version
of R has to be specified. Please note that per default the path is set
to /usr/bin/R.

The C<$edit_control> flag has to be set if information on the edit
distance of the read alignments should be stored by
L<Bio::ViennaNGS::BamStat> . If set, this information can be
visualized in a subsequent step by L<Bio::ViennaNGS::BamStatSummary>.

L<Bio::ViennaNGS::BamStat> and L<Bio::ViennaNGS::BamStatSummary> are
in principle compatible with any BAM file from any read
aligner. Nevertheless, one has to be aware that some mapping tools
differ in BAM dialect or with respect to the information stored in the
BAM file. A special flag to use the auxiliary information stored in a
BAM file produced by I<segemehl> is available as
C<$segemehl_control>. Set to '1' if your input file was mapped with
I<segemehl>, as it the case for the provided C1R1.bam. Otherwise set
to '0'.

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


=head2 Creating a new BamStatSummary object.

 $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						   outpath        => $odir,
						   is_segemehl    => $segemehl_control,
						   control_edit   => $edit_control,
						  );


  Initializes a new BamStatSummary object representing data from
  I<segemehl> mapped BAM files (C<<< is_segemehl => 1 >>>)
  specified in @bams. The output directory is set to $odir ('./'),
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

In the next step the initialized data object is populated with
essential data extracted from each read BAM file.

 $bamsummary->populate_data();

For each BAM file in @bams the method C<< new >> from
L<BIO::ViennaNGS::BamStat> is called like this, 

 $bo = Bio::ViennaNGS::BamStat->new(bam => $bamfile);

=cut

$bamsummary->populate_data();

=head2 Quantify data from $bamsummary

A basic measure in the quantification of BAM files is how many reads
are uniquely or multi mapped and how many alignments exist in total.
Depending on whether single-end or paired-end reads are analyzed, the
number of mapped pairs is essential for the latter. To this end we
compile quantitative information for all reads stored in $bamsummary
by,

 $bamsummary->populate_countStat();

=cut

$bamsummary->populate_countStat();

=head2 Produce output for read quantification.

Next, statistical output will be written to a CSV file in $odir, which
can easily be screened with any text editor or spreadsheet program.

 $bamsummary->dump_countStat("csv");

=cut

$bamsummary->dump_countStat("csv");

=head2 Plot read quantification.

To get a visual overview of the consistency of examined samples,
bar-plots can be produced via,

 $bamsummary->make_BarPlot();

which creates a bar-plot in pdf format of read quantification in $odir.

=cut

$bamsummary->make_BarPlot();

=head2 Plot edit distance distribution.

To gain a quick overview of the quality of different mapped RNA-seq
samples, we plot the distribution of edit distances for all reads
aligned to the reference genome for all samples in @bams.

 $bamsummary->make_BoxPlot("data_edit") if( $bamsummary->has_control_edit );

=cut

$bamsummary->make_BoxPlot("data_edit")  if( $bamsummary->control_edit  && $bamsummary->has_control_edit );

=head2 Summary

We used L<Bio::ViennaNGS::BamStat> and
L<Bio:ViennaNGS::BamStatSummary> to extract, store, summarize, and
visualize quantitative and qualitative data stored in a BAM file. Only
exemplary features of the library were illustrated. Further useful
functions are implemented in the corresponding script
bam_quality_stat.pl, or can be implemented according to ones needs.

=head1 AUTHOR

Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=cut
