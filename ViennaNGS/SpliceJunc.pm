# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-19 00:26:26 mtw>
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael Thomas Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# * This library is free software; you can redistribute it and/or modify
# * it under the same terms as Perl itself, either Perl version 5.12.4 or,
# * at your option, any later version of Perl 5 you may have available.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# *
# ***********************************************************************

package ViennaNGS::SpliceJunc;

use Exporter;
use version; our $VERSION = qv('0.01');
use strict;
use warnings;
use File::Spec; # TODO: perform file name operations with File::Spec
use Data::Dumper;

our @ISA = qw(Exporter);

our @EXPORT = qw(bed6_ss_from_bed12 bed6_ss_from_segemehl_splits
		 novel_sj_from_intersect_annot_segesplit );

our @EXPORT_OK = qw( );


# bed6_ss_from_bed12( $bed12,$dest_dir,$window )
#
# Extracts splice junctions from bed12.
#
# Writes a bed6 file for each transcript found in the bed12, listing
# all splice sites of this transcript, optionally flanking it with a
# window of +/-$window nt.
sub bed6_ss_from_bed12{
  my ($bed12,$dest_dir,$window) = @_;
  my ($pos5,$pos3);
  my $splicesites = 0;
  my @bedline = ();

  die("ERROR [ViennaNGS::SpliceJunc::bed6_ss_from_bed12()] $bed12 does
  not exists\n") unless (-e $bed12);

  die("ERROR [ViennaNGS::SpliceJunc::bed6_ss_from_bed12()] $dest_dir
  does not exist") unless (-d $dest_dir);

  unless ($dest_dir =~ /\/$/){$dest_dir .= "/";}

  #print ">> bed12 -> $bed12\n";
  #print ">> dest_dir -> $dest_dir\n";
  open(BED12IN, "< $bed12") or die $!;
  while(<BED12IN>){
    chomp;
    my ($chr,$chromStart,$chromEnd,$name,$score,$strand,$thickStart,$thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split("\t");
    my @blockSize  = split(/,/,$blockSizes);
    my @blockStart = split(/,/,$blockStarts);
    unless (scalar @blockSize == scalar @blockStart){
      die "ERROR: unequal element count in blockStarts and blockSizes\n";
    }
    my $fn = sprintf("%s_%d-%d_%s.annotatedSS.bed6",$chr,$chromStart,$chromEnd,$name);
    my $bed6_fn = $dest_dir.$fn;
    my $tr_count = 1;

    if ($blockCount >1){ # only transcripts with 2 or more exons (!!!)

      open(BED6OUT, "> $bed6_fn") or die
      "[ViennaNGS::SpliceJunc::bed6_ss_from_bed12() Cannot open
      BED6OUT $!";

      for (my $i=0;$i<$blockCount-1;$i++){
	$pos5 = $chromStart+$blockStart[$i]+$blockSize[$i];
	$pos3 = $chromStart+$blockStart[$i+1];
	#$seqPOS = get_stranded_subsequence($fastaobj{$chr},$pos5+1,$pos5+$ss_motif_length,"+");
	#$seqNEG = get_stranded_subsequence($fastaobj{$chr},$pos3-$ss_motif_length+1,$pos3,"+");
	$splicesites++;
	#if (($seqPOS eq "GT" && $seqNEG eq "AG") ||
	#	  ($seqPOS eq "CT" && $seqNEG eq "AC")) {
	# $cansplicesites++;
	#	#$cansplits += $reads;
	#	$is_canonical="c";
	#}
	#$splicemotif = join("|",$seqPOS,$seqNEG);
	my $tr_name = sprintf("%s.%02d",$name,$tr_count++);
	@bedline = join("\t",$chr,eval($pos5-$window),eval($pos5+$window),$tr_name,$score,$strand);
	print BED6OUT "@bedline\n";
	@bedline = join("\t",$chr,eval($pos3-$window),eval($pos3+$window),$tr_name,$score,$strand);
	print BED6OUT "@bedline\n";
	#print "i=$i\t$pos5\t$pos3\t($seqPOS|$seqNEG)\n";
	#print "@bedline\n";
      } # end for
      close(BED6OUT);
    } # end if
  }
  close(BED12IN);
}

# bed6_ss_from_segemehl_splits ( $segesplitbed,$dest_dir,$window,$mincov )
#
# Extracts splice junctions from segemehl/testrealign -n bed6
#
# Writes a bed6 file for each splice junction found in the
# segemehl/testrealign -n input, optionally flanking it with a window
# of +/-$window nt. Only splice junctions supported by at least $mcov
# reads are considered.
sub bed6_ss_from_segemehl_splits{
  my ($segesplitbed,$dest_dir,$window,$mcov) = @_;
  my ($reads,$proper,$passed);
  my $too_low_coverage = 0;
  my @bedline = ();

  die("ERROR [ViennaNGS::SpliceJunc::bed6_ss_from_segemehl_splits()]
  $segesplitbed does not exist\n") unless (-e $segesplitbed);

  die("ERROR [ViennaNGS::SpliceJunc::bed6_ss_from_segemehl_splits()]
  $dest_dir does not exist") unless (-d $dest_dir);

  unless ($dest_dir =~ /\/$/){$dest_dir .= "/";}

  open(SEGESPLITS, "< $segesplitbed") or die $!;
  while(<SEGESPLITS>){
    chomp;
    my ($chr, $start, $end, $info, $score, $strand) = split("\t");

    if ($info =~ /^splits\:(\d+)\:(\d+)\:(\d+):(\w):(\w)/){
      $reads = $1;
      $proper = $4;
      $passed = $5;
      next unless ($proper eq 'N');
      next unless ($passed eq 'P');
      if($reads < $mcov){
	$too_low_coverage++;
	next;
      }
    }
    else {
      die "unsupported INFO filed in input BED:\n$info\n";
    }

    my $fn = sprintf("%s_%d-%d.segemehlSS.bed6",$chr,$start,$end);
    my $bed6_fn =  $dest_dir.$fn;
    open(BED6OUT, "> $bed6_fn");
    @bedline = join("\t",$chr,eval($start-$window),eval($start+$window),$info,$score,$strand);
    print BED6OUT "@bedline\n";
    @bedline = join("\t",$chr,eval($end-$window),eval($end+$window),$info,$score,$strand);
    print BED6OUT "@bedline\n";
    close(BED6OUT);
  }
  close(SEGESPLITS);
}

# novel_sj_from_intersect_annot_segesplit ( $p_annot,$p_segesplit,$p_out )
#
# Intersect splice junctions determined by segemehl with annotated
# splice junctions. Determine novel and existing splice junctions.
#
# Writes bed6 files for existing and novel splice junctions to $p_out.
sub novel_sj_from_intersect_annot_segesplit{
  my ($p_annot,$p_segesplit,$p_out,$prefix) = @_;
  my ($processed_segesplit_junctions) = 0x1;
  my @segesplit_junctions = ();
  my @transcript_beds = ();
  my %asj = (); # annotated splice junctions hash
  my $bedtools = `which bedtools`; chomp($bedtools);
  my $sortBed  = `which sortBed`; chomp($sortBed);

  die("ERROR
  [ViennaNGS::SpliceJunc::novel_sj_from_intersect_annot_segesplit()]
  $p_annot does not exist\n") unless (-d $p_annot);
  unless ($p_annot =~ /\/$/){$p_annot .= "/";}

  die("ERROR
  [ViennaNGS::SpliceJunc::novel_sj_from_intersect_annot_segesplit()]
  $p_segesplit does not exist\n") unless (-d $p_segesplit);
  unless ($p_segesplit =~ /\/$/){$p_segesplit .= "/";}

  die("ERROR
  [ViennaNGS::SpliceJunc::novel_sj_from_intersect_annot_segesplit()]
  $p_out does not exist\n") unless (-d $p_out);
  unless ($p_out =~ /\/$/){$p_out .= "/";}

  # get a list of all files in $p_annot
  opendir(my $dh_annotated_tr, $p_annot)  or die "Can't opendir $p_annot: $!";
  while(readdir $dh_annotated_tr) {
    push @transcript_beds, $_;
  }
  closedir($dh_annotated_tr);
  #print Dumper(\@transcript_beds);

  # get list of all segemehl splice junctions
  opendir(my $dh_segesplit, $p_segesplit) or die "Can't opendir $p_segesplit: $!";
  @segesplit_junctions = grep { /^(chr\d+)\_(\d+)-(\d+)\.segemehlSS\.bed6/ } readdir($dh_segesplit);

  # process segemehl splice junctions files
  foreach my $file (@segesplit_junctions){
    $processed_segesplit_junctions++;
   # print "processing $file ... ";
    die "Unexpected file name pattern\n" unless ($file =~ /(chr\d+)\_(\d+)-(\d+)/);
    my $sc = $1;
    my $s5 = $2;
    my $s3 = $3;
    my $pattern = $sc."_".$s5."-".$s3;

    my @annotated_beds = grep { /^(chr\d+)\_(\d+)-(\d+)/ && $2<=$s5 && $3>=$s3 && $sc eq $1} @transcript_beds;
    #print"\t intersecting against ".(eval $#annotated_beds+1)." transcripts: @annotated_beds \n";

    # intersect current segemehl SJ against all transcripts in @annotated_beds
    foreach my $i (@annotated_beds){
      my $a = $p_segesplit.$file;
      my $b = $p_annot.$i;
      my $intersect_cmd = "$bedtools intersect -a $a -b $b -c";
      open(INTERSECT, $intersect_cmd."|");
      while(<INTERSECT>){
	chomp;
#	print ">> $_";
	if ($_ =~ /1$/) {
	  $asj{$pattern} = 1;
#	  print " annotated SS";
	}
#	print "\n";
      }
      close(INTERSECT);
    }
    if ( $processed_segesplit_junctions % 1000 == 0){
      print "processed $processed_segesplit_junctions segesplit junctions\n";
    }
  }
  closedir($dh_segesplit);

  # go through the segemehl splice junctions files once more and discriminate
  # novel from existing splice junctions
  if (length($prefix)>0){$prefix .= ".";}
  my $outname_exist = $p_out.$prefix."exist.SS.bed";
  my $outname_exist_u = $p_out.$prefix."exist.SS.u.bed";
  my $outname_novel =  $p_out.$prefix."novel.SS.bed";
  my $outname_novel_u =  $p_out.$prefix."novel.SS.u.bed";
  open (EXISTOUT, "> $outname_exist_u");
  open (NOVELOUT, "> $outname_novel_u");

  # write new ones to NOVELOUT; existing ones to EXISTOUT
  foreach my $file (@segesplit_junctions){
    if ($file =~ m/^(chr\d+\_\d+-\d+)/){
       my $pattern = $1;
       my $fn = $p_segesplit.$file;
       open (SJ, "< $fn") or die "Cannot open $fn $!";
       if (exists $asj{$pattern}){ # annotated splice junction
	 while(<SJ>){ print EXISTOUT "$_";}
       }
       else { # novel splice junction
	 while(<SJ>){ print NOVELOUT "$_";}
       }
       close(SJ);
     }
    else{ warn "Error with parsing segesplit junction file names\n";}
  }
  close(EXISTOUT);
  close(NOVELOUT);

  # sort the resulting bed files
  my $cmdl = "$sortBed -i  $outname_exist_u > $outname_exist";
  system($cmdl);
  unlink($outname_exist_u);
  $cmdl = "$sortBed -i  $outname_novel_u > $outname_novel";
  system($cmdl);
  unlink($outname_novel_u);

  # TODO: Test for canonical / non-canonical motif (write routine for that)
  printf STDERR "processed $processed_segesplit_junctions segemehl junctions\n";
}

1;
__END__

=head1 NAME

ViennaNGS::SpliceJunc - Perl extension for Alternative Splicing Analysis

=head1 SYNOPSIS

  use ViennaNGS::SpliceJunc;

  bed6_ss_from_bed1($bed12,$dest_dir,$window)
  bed6_ss_from_segemehl_splits($segesplitbed,$dest_dir,$window,$mcov)
  novel_sj_from_intersect_annot_segesplit($p_annot,$p_segesplit,$p_out,$prefix)

=head1 DESCRIPTION

ViennaNGS::SpliceJunc is a Perl module for Alternative Splicing (AS)
analysis. It provides routines for identification and characterization of
novel and existing (annotated) splice junctions from RNA-seq data.

Identification of novel splice junctions is based on insecting potentially
novel splice junctions from RNA-seq data with annotated splice junctions.

=head2 EXPORT

Routines:
  bed6_ss_from_bed1($bed12,$dest_dir,$window)
  bed6_ss_from_segemehl_splits(($segesplitbed,$dest_dir,$window,$mcov)
  novel_sj_from_intersect_annot_segesplit($p_annot,$p_segesplit,$p_out,$prefix)

Variables:
   none

=head3 bed6_ss_from_bed1($bed12,$dest_dir,$window)

Extracts splice junctions from an BED12 file (provided via argument
$bed12), writes a BED6 file for each transcript to $dest_dir,
containing all its splice junctions. Output splice junctions can be
flanked by a window of +/- $window nt. Each splice junction is
represented as two bed lines in the output BED6.

=head3 bed6_ss_from_segemehl_splits($segesplitbed,$dest_dir,$window,$mcov)

Extracts splice junctions from segemehl's haarz / testrealign -n BED6
output and writes a BED6 file for each splice junction given in the
input to $dest_dir. Output splice junctions can be flanked by a window
of +/- $window nt. Each splice junction is represented as two bed
lines in the output BED6.

=head3 novel_sj_from_intersect_annot_segesplit($p_annot,$p_segesplit,$p_out,$prefix)

Intersects all splice junctions identified in an RNA-seq experiment
with annotated splice junctions. Identifies and characterized novel
and existing splice junctions. Each BED6 file in $p_segesplit is
intersected with those transcript splice junction BED6 files in
$P_annot, whose genomic location spans the query splice junction. This
is just to prevent the tool from intersecting each splice site found
in the mapping data with all annotated transcripts.

The intersection operations are performed with intersectBed from the
BEDtools suite (https://github.com/arq5x/bedtools2). BED sorting
operations are performed with sortBed.

Writes two BEd6 files to $p_out (optionally prefixed by $prefix),
which contain novel and existing splice junctions respectively.

=head1 DEPENDENCIES

ViennaNGS::SpliceJunc uses third-party tools for computing intersections of
BED files. Specifically, the intersectBed utility from the BEDtools suite
(http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html) is
used to compute overlaps and sortBed is used to sort BED output files
(http://bedtools.readthedocs.org/en/latest/content/tools/sort.html). Make
sure that those third-party utilities are available on your system, and
that hey can be found by the perl interpreter. We recommend installing the
latest version of BEDtools (https://github.com/arq5x/bedtools2) on your
system.


=head1 SEE ALSO

  perldoc ViennaNGS
  perldoc ViennaNGS::AnnoC

=head1 AUTHOR

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Michael Thomas Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.


=cut
