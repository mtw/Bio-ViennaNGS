# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-25 17:02:12 mtw>
#
# TODO:
#       - extractfilter canonical splice junctions in intersect_sj()
#       - create bigBed output for novel/exist splice junctions

package ViennaNGS::SpliceJunc;

use Exporter;
use version; our $VERSION = qv('0.02_01');
use strict;
use warnings;
use File::Spec; # TODO: perform file name operations with File::Spec
use Data::Dumper;
use ViennaNGS;

our @ISA = qw(Exporter);

our @EXPORT = qw(bed6_ss_from_bed12 bed6_ss_from_rnaseq
		 intersect_sj );

our @EXPORT_OK = qw( );


# bed6_ss_from_bed12( $bed12,$dest_dir,$window,$fastaobjR )
#
# Extracts splice junctions from bed12.
#
# Writes a BED6 file for each transcript found in the BED12, listing
# all splice sites of this transcript, optionally flanking it with a
# window of +/-$window nt.
#
sub bed6_ss_from_bed12{
  my ($bed12,$dest_dir,$window,$fastaobjR) = @_;
  my ($i,$c,$tr_name,$pos5,$pos3);
  my $splicesites = 0;
  my @bedline = ();
  my $this_function = (caller(0))[3];

  die("ERROR [$this_function] $bed12 does not exists\n")
    unless (-e $bed12);

  die("ERROR [$this_function] $dest_dir does not exist")
    unless (-d $dest_dir);

  unless ($dest_dir =~ /\/$/){$dest_dir .= "/";}

  open(BED12IN, "< $bed12") or die $!;
  while(<BED12IN>){
    chomp;
    my ($chr,$chromStart,$chromEnd,$name,$score,$strand,$thickStart,
	$thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split("\t");
    my @blockSize  = split(/,/,$blockSizes);
    my @blockStart = split(/,/,$blockStarts);
    unless (scalar @blockSize == scalar @blockStart){
      die "ERROR: unequal element count in blockStarts and blockSizes\n";
    }
    my $fn = sprintf("%s_%d-%d_%s.annotatedSS.bed6",$chr,$chromStart,$chromEnd,$name);
    my $bed6_fn = $dest_dir.$fn;
    my $tr_count = 1;

    if ($blockCount >1){ # only transcripts with 2 or more exons (!!!)

      open(BED6OUT, "> $bed6_fn") or die "cannot open BED6OUT $!";

      for ($i=0;$i<$blockCount-1;$i++){
	$pos5 = $chromStart+$blockStart[$i]+$blockSize[$i];
	$pos3 = $chromStart+$blockStart[$i+1];
	$c = ss_isCanonical($chr,$pos5,$pos3,$fastaobjR);
	$tr_name = sprintf("%s.%02d",$name,$tr_count++);
	@bedline = join("\t",$chr,eval($pos5-$window),
			eval($pos5+$window),$tr_name,$c,$strand);
	print BED6OUT "@bedline\n";
	@bedline = join("\t",$chr,eval($pos3-$window),
			eval($pos3+$window),$tr_name,$c,$strand);
	print BED6OUT "@bedline\n";
      } # end for

      close(BED6OUT);
    } # end if
  } # end while
  close(BED12IN);
}

# bed6_ss_from_rnaseq ( $inbed,$dest_dir,$window,$mincov )
#
# Extracts splice junctions from mapped RNA-seq data
#
# Writes a BED6 file for each splice junction present in the input,
# optionally flanking it with a window of +/-$window nt. Only splice
# junctions supported by at least $mcov reads are considered.
sub bed6_ss_from_rnaseq{
  my ($inbed,$dest_dir,$window,$mcov,$fastaobjR) = @_;
  my ($reads,$proper,$passed,$c,$pos5,$pos3);
  my $too_low_coverage = 0;
  my @bedline = ();
  my $this_function = (caller(0))[3];

  die("ERROR [$this_function] $inbed does not exist\n")
    unless (-e $inbed);

  die("ERROR [$this_function] $dest_dir does not exist")
    unless (-d $dest_dir);

  unless ($dest_dir =~ /\/$/){$dest_dir .= "/";}

  open(INBED, "< $inbed") or die $!;
  while(<INBED>){
    chomp;
    my ($chr, $start, $end, $info, $score, $strand) = split("\t");
    $end = $end-2; # required for segemehl/testrealign BED6 files

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
    $pos5 = $start;
    $pos3 = $end;
    $c = ss_isCanonical($chr,$pos5,$pos3,$fastaobjR);

    my $fn = sprintf("%s_%d-%d.mappedSS.bed6",$chr,$start,$end);
    my $bed6_fn =  $dest_dir.$fn;
    open(BED6OUT, "> $bed6_fn");
    @bedline = join("\t",$chr,eval($start-$window),
		    eval($start+$window),$info,$c,$strand);
    print BED6OUT "@bedline\n";
    @bedline = join("\t",$chr,eval($end-$window),
		    eval($end+$window),$info,$c,$strand);
    print BED6OUT "@bedline\n";
    close(BED6OUT);
  }
  close(INBED);
}

# intersect_sj($p_annot,$p_mapped,$p_out,$prefix,$window,$mil,$can)
#
# Intersect splice junctions determined by RNA-seq with annotated
# splice junctions. Determine novel and existing splice junctions.
#
# Writes BED6 files for existing and novel splice junctions to $p_out.
sub intersect_sj{
  my ($p_annot,$p_mapped,$p_out,$prefix,$window,$mil,$can) = @_;
  my ($dha,$dhm);
  my $processed = 0;
  my @junctions = ();
  my @transcript_beds = ();
  my %asj = (); # annotated splice junctions hash
  my $bedtools = `which bedtools`; chomp($bedtools);
  my $sortBed  = `which sortBed`; chomp($sortBed);
  my $this_function = (caller(0))[3];

  die("ERROR [$this_function] $p_annot does not exist\n") unless (-d $p_annot);
  unless ($p_annot =~ /\/$/){$p_annot .= "/";}

  die("ERROR [$this_function] $p_mapped does not exist\n") unless (-d
  $p_mapped);
  unless ($p_mapped =~ /\/$/){$p_mapped .= "/";}

  die("ERROR [$this_function] $p_out does not exist\n") unless (-d
  $p_out);
  unless ($p_out =~ /\/$/){$p_out .= "/";}

  # get a list of all files in $p_annot
  opendir($dha, $p_annot)  or die "Can't opendir $p_annot: $!";
  while(readdir $dha) { push @transcript_beds, $_; }
  closedir($dha);

  # get a list of all splice junctions seen in RNA-seq data
  opendir($dhm, $p_mapped) or die "Can't opendir $p_mapped: $!";
  @junctions = grep { /^(chr\d+)\_(\d+)-(\d+)\.mappedSS\.bed6/ } readdir($dhm);

  # process splice junctions seen in RNA-seq data
  foreach my $file (@junctions){
    $processed++;
    # print "processing $file ... ";
    die "Unexpected file name pattern\n" unless ($file =~ /(chr\d+)\_(\d+)-(\d+)/);
    my $sc = $1;
    my $s5 = $2;
    my $s3 = $3;
    my $pattern = $sc."_".$s5."-".$s3;

    my @annotated_beds = grep { /^(chr\d+)\_(\d+)-(\d+)/ && $2<=$s5 && $3>=$s3 && $sc eq $1} @transcript_beds;
    #print"\t intersecting against ".(eval $#annotated_beds+1)." transcripts: @annotated_beds \n";

    # intersect currently opened SJ against all transcripts in @annotated_beds
    foreach my $i (@annotated_beds){
      my $a = $p_mapped.$file;
      my $b = $p_annot.$i;
      my $intersect_cmd = "$bedtools intersect -a $a -b $b -c -nobuf";
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
    if ( $processed % 1000 == 0){
      print "processed $processed splice junctions\n";
    }
  }
  closedir($dhm);

  # go through the mapped splice junctions files once more and
  # separate novel from existing splice junctions
  if (length($prefix)>0){$prefix .= ".";}
  my $outname_exist = $p_out.$prefix."exist.SS.bed";
  my $outname_exist_u = $p_out.$prefix."exist.SS.u.bed";
  my $outname_novel =  $p_out.$prefix."novel.SS.bed";
  my $outname_novel_u =  $p_out.$prefix."novel.SS.u.bed";
  open (EXISTOUT, "> $outname_exist_u");
  open (NOVELOUT, "> $outname_novel_u");

  # write new ones to NOVELOUT; existing ones to EXISTOUT
  foreach my $file (@junctions){
    if ($file =~ m/^(chr\d+\_\d+-\d+)/){
       my $pattern = $1;
       my $fn = $p_mapped.$file;
       open (SJ, "< $fn") or die "Cannot open $fn $!";
       while(<SJ>){
	 chomp;
	 $_ = m/^(chr\w+)\s(\d+)\s(\d+)\s(splits:\d+:\d+:\d+:\w:\w)\s(\d+)\s([+-01])/;
	 my $chr = $1;
	 my $start = $2;
	 my $end = $3;
	 my $name = $4;
	 my $score = $5;
	 my $strand = $6;
	 my @bedline = join("\t",$chr,eval($start+$window),eval($start+$window+1),$name,$score,$strand);
	 if (exists $asj{$pattern}){ # annotated splice junction
	   print EXISTOUT "@bedline\n";
	 }
	 else {  # novel splice junction
	   print NOVELOUT "@bedline\n";
	 }
       } # end while
       close(SJ);
     } # end if
    else{ warn "Error with parsing BED6 junction file names __FILE__ __LINE__\n";}
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
  printf STDERR "processed $processed splice junctions\n";
}

# ss_isCanonical ( $chr,$p5,$p3,$fastaobjR )
#
# Possible scenarios for canonical splice junction (option 1 shown
# here)
#   ------------------>
# 5'===]GU..........AG[====3'
#
#   <-----------------
# 3'===]GA..........UG[====5'
#
# in terms of the underlying reference sequence this boils down to the following combinations
# 5'===]GT|CT....AG|AC[====3' ie GT->AG or CT->AC
# 5'===]GC|CT....AG|GC[====3' ie GC->AG or CT->GC
# 5'===]AT|GT....AC|AT[====3' ie AT->AC or GT->AT
sub ss_isCanonical{
  my ($chr,$p5,$p3,$fastaobjR) = @_;
  my ($seqL,$seqR,$pattern);
  my %fastaO = %$fastaobjR;
  my $ss_motif_length = 2;
  my $this_function = (caller(0))[3];
  my $c = -1;

  # am besten hier die nt Kombinatioinen irgendwie definieren, damit wir die length dann unten haben
  $seqL = get_stranded_subsequence($fastaO{$chr},$p5+1,$p5+$ss_motif_length,"+");
  $seqR = get_stranded_subsequence($fastaO{$chr},$p3-$ss_motif_length+1,$p3,"+");

  $pattern = sprintf("%s|%s",$seqL,$seqR);
  #print STDERR "[$this_function] p5->p3 ($p5 -- $p3) $seqL -- $seqR\n";

  if ($pattern =~ /^GT|AG$/) { $c = 1000;}
  elsif ($pattern =~ /^CT|AC$/) { $c = 1000; }
  elsif ($pattern =~ /^GC|AG$/) { $c = 1000; }
  elsif ($pattern =~ /^CT|GC$/) { $c = 1000; }
  elsif ($pattern =~ /^AT|AC$/) { $c = 1000; }
  elsif ($pattern =~ /^GT|AT$/) { $c = 1000; }
  else { $c = 0;}

  return $c;
}

1;
__END__

=head1 NAME

ViennaNGS::SpliceJunc - Perl extension for alternative splicing analysis

=head1 SYNOPSIS

  use ViennaNGS::SpliceJunc;

  bed6_ss_from_bed12($bed12,$dest_dir,$window,$fastaobjR)
  bed6_ss_from_rnaseq($inbed,$dest_dir,$window,$mcov)
  intersect_sj($p_annot,$p_mapped,$p_out,$prefix,$window,$mil,$can)

=head1 DESCRIPTION

ViennaNGS::SpliceJunc is a Perl module for alternative splicing (AS)
analysis. It provides routines for identification and characterization of
novel and existing (annotated) splice junctions from RNA-seq data.

Identification of novel splice junctions is based on insecting potentially
novel splice junctions from RNA-seq data with annotated splice junctions.

=head2 EXPORT

Routines:
  bed6_ss_from_bed12($bed12,$dest_dir,$window,$fastaobjR)
  bed6_ss_from_rnaseq($inbed,$dest_dir,$window,$mcov)
  intersect_sj($p_annot,$p_mapped,$p_out,$prefix,$window,$mil,$can)

Variables:
   none

=head1 SUBROUTINES

=head2 bed6_ss_from_bed12($bed12,$dest_dir,$window,$fastaobjR)

Extracts splice junctions from an BED12 file (provided via argument
$bed12), writes a BED6 file for each transcript to $dest_dir,
containing all its splice junctions. Output splice junctions can be
flanked by a window of +/- $window nt. $fastaobjR is a reference to a
Bio::PrimarySeq::Fasta object holding the underlying reference
genome. Each splice junction is represented as two bed lines in the
output BED6.

=head2 bed6_ss_from_rnaseq($inbed,$dest_dir,$window,$mcov)

Extracts splice junctions from mapped RNA-seq data. The input BED6
file should contain coordinates of introns in the following syntax:

chr1    3913    3996    splits:97:97:97:N:P     0       +

The fourth column in this BED file (correponding to the 'name' field)
should be a colon-separated string of six elements, where the first
element should be 'splits' and the second element is assumed to hold
the number of reads supporting this splice junction. The fifth element
indicates the splice junction type: A capital 'N' determines a normal
splice junction, whereas 'C' indicates circular and 'T' indicates
trans-splice junctions, respectively. Only normal splice junctions
('N') are considered, the rest is skipped. Elements 3, 4 and 6 are not
further processed.

We recommend using I<segemehl> for generating this type of BED6
files. This routine is, however, not limited to segemehl output. BED6
files containing splice junction information from other short read
mappers or third-party sources will be processed if hey are formatted
as described above.

This routine writes a BED6 file for each splice junction provided in
the input to $dest_dir. Output splice junctions can be flanked by a
window of +/- $window nt. Each splice junction is represented as two
bed lines in the output BED6.

=head2 intersect_sj($p_annot,$p_mapped,$p_out,$prefix,$window,$mil,$can)

Intersects all splice junctions identified in an RNA-seq experiment
with annotated splice junctions. Identifies and characterized novel
and existing splice junctions. Each BED6 file in $p_mapped is
intersected with those transcript splice junction BED6 files in
$p_annot, whose genomic location spans the query splice junction. This
is just to prevent the tool from intersecting each splice site found
in the mapped RNA-seq data with all annotated transcripts. $mil
specifies a maximum intron length. $can is a switch for filtering
canonical splice junctions.

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

=head1 AUTHORS

Michael Thomas Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.


=cut
