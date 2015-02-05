#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-01-27 16:24:20 fabian>

=head1 NAME

Tutorial_pipeline00.pl - Infering qualkity parameter form a BAM file.

=head1 SYNOPSIS

  Tutorial_pipeline00.pl

=head1 DESCRIPTION

Tutorial how to use Bio::ViennaNGS::BamStat and
Bio::ViennaNGS::BamStatSummary to retrieve quality and quantitiy
statistics from input BAM file.

=head1 Prcedure

=head2 Include libraries

 use Data::Dumper;
 use Bio::ViennaNGS::BamStat;
 use Bio::ViennaNGS::BamStatSummary;

=over 4

=item * Data::Dumper

For easy access to complex data structure.

=item * Bio::ViennaNGS::BamStat

Extracts quality and quantity paramters from an input BAM file.

=item * Bio::ViennaNGS::BamStatSummary

Sumarizes, compares and plots data compiled by Bio::ViennaNGS::BamStat.

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

=back

=head2 Define control variables

 @bams     = qw# C1R2.bam #;
 $odir     = '.';
 $rlibpath = '/usr/bin/R';

 $edit_control     = 1;
 $segemehl_control = 1;


=over 4

=item * @bams

Array with all BAM files, including their path, to be processed. In
the course of this tutorial please retrieve the C1R1.bam file from
xxxxxx and store it in the current working directory. Other inputs
will not be accepted.

=item * $odir

Path to the directory where the output files will be created. If files
with same names do already exist in this particular directory, they
will be overwritten.

=item * $rlibpath

Path to the installation of R.

=item * $edit_control

Control flag. Set to 1 if a statistics of the edit distance of each
read should be reported. Otherwise set to 0;

=item * $segemehl_control

Control flag. Set to 1 if the input BAM file was produced by the short
read mapper I<segemehl>. Takes care of segemehl specific BAM dialect
issues. Otherwise set to 0;

=cut

my @bams     = qw# /home/mescalin/fabian/Work/ViennaNGS/Data/paired-end/segemehl_test1.bam #;
#next unless ($bam[0] eq $bam[-1] && $bam[0] eq "C1R1.bam");
my $odir     = '/home/mescalin/fabian/Work/ViennaNGS/Tutorial/Prog';
my $rlibpath = '/usr/bin/R';
my %data     = ();

=back

=head2 Creating new BamStatSummary object.

 $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						   outpath        => $odir,
						   rlib           => $rlibpath,
						   is_segemehl    => 1,
						   control_edit   => 1,
						  );

=over 4

=item * Options

Initialize new BamStatSummary object representing data from all
segemehl BAM files in @bams, setting the output directory to $odir,
where beside standard read quantification also the edit distance of
each read will be reported.

=cut


my $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						     outpath        => $odir,
						     rlib           => $rlibpath,
                                                     is_segemehl    => 1,
                                                     control_edit   => 1,
                                                     control_flag   => 1,
                                                     control_score  => 1,

						    );


=back

=head2 Read-in BAM files @bams

Processes each BAM file in @bams and compiles the relevant data into $bamsummary.

 $bamsummary->populate_data();

Thereby, for each BAM file in @bams the method new from BIO::ViennaNGS::BamStat is called like this

 $bo = Bio::ViennaNGS::BamStat->new(bam => $bamfile);

Use C<print Dumper($bamsummary);> to check the object.

=cut

$bamsummary->populate_data();

=head2 Quantify data from $bamsummary

Compiles quantitative information for all reads stored in $bamsummary.

 $bamsummary->populate_countStat();

Use C<print Dumper($bamsummary);> to check the object.

=cut

$bamsummary->populate_countStat();

=head2 Produce output for read quantification.

Create file for the read quantification in $odir. File formate is
*.csv.

 $bamsummary->dump_countStat("csv");

=cut

$bamsummary->dump_countStat("csv");

=head2 Plot read quantification.

Create a barplot for the read quantification in $odir. File formate is
*.pdf.

 $bamsummary->make_BarPlot();

=cut

$bamsummary->make_BarPlot();

=head2 Plot edit distance distribution.

Create a boxplot of the distribution of edit distances for all reads
and all samples in @bams.

 $bamsummary->make_BoxPlot("data_edit") if($bamsummary->control_edit  && $bamsummary->has_control_edit);

=cut

$bamsummary->make_BoxPlot("data_edit")  if($bamsummary->control_edit  && $bamsummary->has_control_edit);


=head1 AUTHOR

Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=cut
