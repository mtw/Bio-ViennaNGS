#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-02-12 23:31:38 mtw>

=head1 NAME

Tutorial_pipeline00.pl - Inferring quantitative and qualitative
parameters from an input BAM file.

=head1 SYNOPSIS

  perl Tutorial_pipeline00.pl

=head1 DESCRIPTION

This tutorial demonstrates the functionality of
L<Bio::ViennaNGS::BamStat> and L<Bio::ViennaNGS::BamStatSummary> for
inferring mapping statistics from BAM files. We will collect the
number of mapped reads, read alignments, single-end and paired-end
reads and check/compare quality features stored in a BAM file. The
latter is exemplified by using the module to deduce the distribution
of edit distances for all read alignments.

We are going to examine the amount of different read types in an exact
manner and visualize the distribution of edit distances between
aligned reads and the reference genome.

=head2 PREREQUISITES

For running this tutorial on your machine you will need a full
installation of the L<Bio::ViennaNGS> distribution, including all
third party dependencies, i.e. a recent version of
L<bedtools|https://github.com/arq5x/bedtools2> as well as the
F<bedGraphToBigWig> utility from the UCSC source distribution
(available L<here|http://hgdownload.cse.ucsc.edu/admin/exe/>). In
addition, the following input files (which can be downloaded
L<here|http://rna.tbi.univie.ac.at/ViennaNGS>) will be used throughout
this tutorial:

=over

=item F<C1R1.bam>

=back

We will assume that all input files are available and accessible in
your current working directory.

=head2 DISCLAIMER

This tutorial does not cover all features of
L<Bio::ViennaNGS::BamStat> and L<Bio::ViennaNGS::BamStatSummary>. It
is merely meant to illustrate the basic usage principles. For more
details we refer to the documentation of L<Bio::ViennaNGS::BamStat>
and L<Bio::ViennaNGS::BamStatSummary>.

This tutorial works on a real-world biological data set of several
gigabytes in size i.e. the analysis will eventually take a few hours
to finish, depending on your hardware. If you run this script locally
you need to ensure that your system has enough hardware resources
available.

=head1 PROCEDURE

=head2 Include libraries

 use Bio::ViennaNGS::BamStat;
 use Bio::ViennaNGS::BamStatSummary;

The following libraries are included: L<Bio::ViennaNGS::BamStat>
provides methods for extracting key aspects in terms of quality and
quantity information from a BAM file and stores essential information
in a data object. L<Bio::ViennaNGS::BamStatSummary> provides methods
to compare and visualize data stored in the L<Bio::ViennaNGS::BamStat>
object.

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
are set. The array C<@bams> holds a list of all BAM files that are
going to be analyzed. This tutorial is restricted to a single file
F<C1R1.bam>, which should be accessible in the current working
directory.

Since the F<Tutorial_pipeline00.pl> script produces several output
files with fixed file names, an output directory has to be
specified. This is done in the C<$odir> variable, which is set to the
current working directory per default. Please not that files with the
same name in this directory will be overwritten.

Some methods in L<Bio::ViennaNGS::BamStatSummary> depend on the
C<Statistics::R> library. Therefore, a valid path to a working version
of the I<R statistics software> has to be specified. Please note that
the default path is set to /usr/bin/R.

The C<$edit_control> flag has to be set if information on the edit
distance of the read alignments should be stored by
L<Bio::ViennaNGS::BamStat>. If set, this information can be
visualized in a subsequent step by L<Bio::ViennaNGS::BamStatSummary>.

L<Bio::ViennaNGS::BamStat> and L<Bio::ViennaNGS::BamStatSummary> are
(in principle) compatible with any BAM file produced by any read
aligner. Nevertheless, one has to be aware that some mappers produce
different BAM dialects that differ in the type and amount of
information stored in the BAM file. A special flag to use the
auxiliary information stored in a BAM file produced by I<segemehl> is
available as C<$segemehl_control>. Set to '1' if your input file has
been produced by I<segemehl>, as it the case for the provided
F<C1R1.bam>. Otherwise set to '0'.

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

 my $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						   outpath        => $odir,
						   is_segemehl    => $segemehl_control,
						   control_edit   => $edit_control,
						  );


Initializes a new L<Bio::ViennaNGS::BamStatSummary> object that
contains data from BAM files maped by I<segemehl> (C<is_segemehl
=E<gt> 1>), specified in @bams. The output directory is set to $odir
('./'), where standard read quantification and the edit distance
(C<control_edit =E<gt> 1>) of each read will be stored.

=cut

my $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						     outpath        => $odir,
                                                     is_segemehl    => 1,
                                                     control_edit   => 1,
                                                     control_flag   => 1,
                                                     control_score  => 1,

						    );


=head2 Read in BAM files @bams

In the next step the initialized data object is populated with
essential data extracted from each read BAM file.

 $bamsummary->populate_data();

For each BAM file in @bams the method C<new> from
L<BIO::ViennaNGS::BamStat> is called like this,

 $bo = Bio::ViennaNGS::BamStat->new(bam => $bamfile);

=cut

$bamsummary->populate_data();

=head2 Quantify data from $bamsummary

A basic measure in the quantification of BAM files is the number of
uniquely or multi mapped reads and how many alignments exist in total.
Depending on whether single-end or paired-end reads are analyzed, the
number of mapped pairs is essential for the latter. To this end we
compile quantitative information for all reads stored in $bamsummary
by

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

To get a visual overview of the consistency of examined samples, bar
plots can be produced via,

 $bamsummary->make_BarPlot();

which creates a barplot in pdf format of read quantification in $odir.

=cut

$bamsummary->make_BarPlot();

=head2 Plot edit distance distribution.

To get a quick overview of the quality of different mapped RNA-seq
samples, we plot the distribution of edit distances for all reads
aligned to the reference genome for all samples in @bams.

 $bamsummary->make_BoxPlot("data_edit") if( $bamsummary->has_control_edit );

=cut

$bamsummary->make_BoxPlot("data_edit")  if( $bamsummary->control_edit  && $bamsummary->has_control_edit );

=head2 Summary

We applied L<Bio::ViennaNGS::BamStat> and
L<Bio:ViennaNGS::BamStatSummary> to extract, store, summarize, and
visualize quantitative and qualitative data stored in a BAM file. Only
exemplary features of the library were illustrated. The script
L<bam_quality_stat.pl|http://search.cpan.org/dist/Bio-ViennaNGS/scripts/bam_quality_stat.pl>.

=head1 AUTHOR

=over

=item Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=item Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=back

=cut
