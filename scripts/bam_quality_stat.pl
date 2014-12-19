#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-19 23:37:35 mtw>

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Path::Class;
use Bio::ViennaNGS::BamStatSummary;


###########################
## # #   Variables   # # ##
###########################

my $bam_dir  = '';           # keeps path to directory with bams
my @bams     = '';           # array holding all bam files to process
my $odir     = '';           # keeps path to dir where results are stored
my $rlibpath = '/usr/bin/R'; # path to R installation (`which R`)
my %data     = ();           # stores all results from single bam files
my $VERBOSE  = 0;

my $match_control    = 0;   # provides stats how many mapped bases match the reference genome
my $clip_control     = 0;   # provides stats how many bases are soft or hard clipped
my $split_control    = 0;   # provides stats how many mapped reads are splitted (and how often);
my $qual_control     = 0;   # provides stats on quality of the match
my $edit_control     = 0;   # provides stats on the edit distance between read and mapped reference
my $flag_control     = 1;   # analyses the sam bit flag for qual/strands/pair_vs_single reads
my %raw_flag_data    = ();
my $score_control    = 1;   # provides stats on per-base quality scores
my $uniq_control     = 1;   # gives number and stats of multiplicity of readaligments (must be =1)
my $segemehl_control = 0;   # toggles to consider segemehl specific bam feature

##################################
## # # Command-line Options # # ##
##################################

pod2usage(-verbose => 0)
        unless GetOptions(
	  "dir|d=s"     => \$bam_dir,
	  "bam|b=s"     =>  sub {@bams=split/:/, $_[1]},
	  "odir|o=s"    => \$odir,
	  "rlib|r=s"    => \$rlibpath,
	  "match!"      => \$match_control,
	  "clip!"       => \$clip_control,
	  "qual!"       => \$qual_control,
	  "edit!"       => \$edit_control,
	  "segemehl!"   => \$segemehl_control,
	  "help|h"      => sub{pod2usage(-verbose => 1)},
	  "man|m"       => sub{pod2usage(-verbose => 2)},
	  "verbose"     => sub{ $VERBOSE++ }
    );

#####################
## # #  Input  # # ##
#####################


if(@bams){
  foreach my $bam (@bams){
    unless (-f $bam){
      warn "Could not find input bam file $bam given via --bam option\n";
      pod2usage(-verbose => 0);

    }
  }
}

if ($bam_dir){
  unless ($bam_dir=~ /^\// || $bam_dir =~ /\.\//){$bam_dir = "./".$bam_dir;}
  unless (-d $bam_dir){
    warn "Could not find input dir $bam_dir given via --dir option\n";
    pod2usage(-verbose => 0);
  }

  opendir(my $dh, $bam_dir) || die "can't opendir $bam_dir: $!\n";
  push @bams,  map {"$bam_dir/$_"} grep { /\.bam$/ && -f "$bam_dir/$_" } readdir($dh);
}

unless ($odir =~ /\/$/){$odir .= "/";}
unless (-d $odir){mkdir $odir or die $!;}

####################
## # #  Main  # # ##
####################

my $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files          => \@bams,
						     outpath        => $odir,
						     rlib           => $rlibpath,
						     is_segemehl    => $segemehl_control,
						     control_match  => $match_control,
						     control_clip   => $clip_control,
						     control_split  =>  $split_control,
						     control_qual   => $qual_control,
						     control_edit   => $edit_control,
						     control_flag   => $flag_control,
						     control_score  => $score_control,
						    );


$bamsummary->populate_data();
$bamsummary->populate_countStat();
$bamsummary->dump_countStat("csv");
$bamsummary->make_BarPlot();


__END__

#####################
## # #   POD   # # ##
#####################


=head1 NAME

bam_quality_stats.pl -- Generates quality statistics from specified
BAM input files.

=head1 SYNOPSIS

bam_quality_stats.pl --dir <PATH> --bam <file1:file2:..> --odir <PATH>
[--match] [--clip] [--qual] [--edit] [--rlib <PATH>] [[--help]

=head1 OPTIONS

=over
    
=item B<--dir>

Path to directory with BAM files. All BAM files in this directory will
be processed.

=item B<--bam>

List seperated by ':'. Specifies bam files which schould be parsed and
processd. Can be combined with --dir.

=item B<--odir>

Path to output directory. In this directory several output files will
be created. Existing files with identiacal names will be over
written. Can be combined with --bam.

=item B<--match>

Provides stats how many bases match the reference genome in the
alignment's CIGAR string. If sequence mismatchs are defined as
alignment matches depends on the used read mapper. Some use the '='
and the 'X' symbol. In this case only sequencing matches are
counted. If the aligner only reports 'M' no distinction between
sequence match and mismatch can be done. In any case, clipped bases,
deletions, insertions are excluded.  =item B<--clip>

Provides stats how many bases are soft or hard clipped.

=item B<--qual>

Provides stats on quality of the read match against the reference
sequence.

=item B<--edit>

Provides stats on the edit distance between read and mapped reference position. 

 =item B<--segemehl>

Toggles to specific segemehl bam file characteristics (e.g., mapped
fragments per line).  =item B<--rlib -r>

Path to the R library. Default is '/usr/bin/R'.

=item B<--help -h>

Print short help.

=item B<--man>

Prints the manual page and exits.

=item B<--verbose>

Prints additional output.

=back

=head1 DESCRIPTION

This program screens a given directory for BAM files and computes
various statistics on them.

=head1 AUTHOR

Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=cut
