# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:33:10 mtw>

package Bio::ViennaNGS::SpliceJunc;

use Exporter;
use version; our $VERSION = qv('0.12_07');
use strict;
use warnings;
use Data::Dumper;
use Bio::ViennaNGS;
use Bio::ViennaNGS::Fasta;
use IPC::Cmd qw(can_run run);
use File::Basename;
use Path::Class;
use Carp;

our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(bed6_ss_from_bed12 bed6_ss_from_rnaseq
		 bed6_ss_to_bed12 intersect_sj ss_isCanonical);

# bed6_ss_from_bed12( $bed12,$dest,$window,$can,$fastaO)
#
# Extracts splice junctions from BED12 annotation.
#
# Writes a BED6 file for each transcript found in the BED12, listing
# all splice sites of this transcript, optionally flanking it with a
# window of +/-$window nt.
sub bed6_ss_from_bed12{
  my ($bed12,$dest,$window,$can,$fastaobjR) = @_;
  my ($i,$tr_name,$pos5,$pos3);
  my ($splicesites,$c,$totalsj,$cansj) = 0x4;
  my @bedline = ();
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] $bed12 does not exists\n"
    unless (-e $bed12);
  croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);

  open(BED12IN, "< $bed12") or croak $!;
  while(<BED12IN>){
    chomp;
    my ($chr,$chromStart,$chromEnd,$name,$score,$strand,$thickStart,
	$thickEnd,$itemRgb,$blockCount,$blockSizes,$blockStarts) = split("\t");
    my @blockSize  = split(/,/,$blockSizes);
    my @blockStart = split(/,/,$blockStarts);
    unless (scalar @blockSize == scalar @blockStart){
      croak "ERROR: unequal element count in blockStarts and blockSizes\n";
    }
    my $fn = sprintf("%s_%d-%d_%s.annotatedSS.bed6",$chr,$chromStart,$chromEnd,$name);
    my $bed6_fn = file($dest,$fn);
    my $tr_count = 1;

    if ($blockCount >1){ # only transcripts with 2 or more exons (!!!)

      open(BED6OUT, "> $bed6_fn") or croak "cannot open BED6OUT $!";

      for ($i=0;$i<$blockCount-1;$i++){
	$totalsj++;
	$pos5 = $chromStart+$blockStart[$i]+$blockSize[$i];
	$pos3 = $chromStart+$blockStart[$i+1];
	if($can){
	  $c = ss_isCanonical($chr,$pos5,$pos3,$fastaobjR);
	  $cansj++;
	}
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

# bed6_ss_from_rnaseq ($bed_in,$dest,$window,$mincov,$can,$fastaO)
#
# Extracts splice junctions from mapped RNA-seq data
#
# Writes a BED6 file for each splice junction present in the input,
# optionally flanking it with a window of +/-$window nt. Only splice
# junctions supported by at least $mcov reads are considered.
sub bed6_ss_from_rnaseq{
  my ($bed_in,$dest,$window,$mcov,$can,$fastaobjR) = @_;
  my ($reads,$proper,$passed,$pos5,$pos3);
  my $c = 0;
  my @bedline = ();
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] $bed_in does not exist\n"
    unless (-e $bed_in);
  croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);

  open(INBED, "< $bed_in") or croak $!;
  while(<INBED>){
    chomp;
    my ($chr, $start, $end, $info, $score, $strand) = split("\t");
    $end = $end-2; # required for segemehl's BED6 files

    if ($info =~ /^splits\:(\d+)\:(\d+)\:(\d+):(\w):(\w)/){
      $reads = $1;
      $proper = $4;
      $passed = $5;
      next unless ($proper eq 'N');
      next unless ($passed =~ /[PFM]$/);
      next if($reads < $mcov);
    }
    else {
      croak "ERROR [$this_function] unsupported INFO field in input BED:\n$info\n";
    }
    $pos5 = $start;
    $pos3 = $end;
    if($can){$c = ss_isCanonical($chr,$pos5,$pos3,$fastaobjR);}
    my $fn = sprintf("%s_%d-%d.mappedSS.bed6",$chr,$start,$end);
    my $bed6_fn = file($dest,$fn);
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

# bed6_ss_to_bed12 ($bed_in,$dest,$window,$mcov)
#
# Produce BED12 from BED6 file holdig splice junctions from mapped
# RNA-seq data
sub bed6_ss_to_bed12{
  my ($bed_in,$dest,$window,$mcov,$circ) = @_;
  my ($reads,$proper,$passed,$pos5,$pos3,$basename,$fn,$bed12_fn);
  my @result = ();
  my @bedline = ();
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] $bed_in does not exist\n"
    unless (-e $bed_in);
  croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);

  $basename = fileparse($bed_in,qr/\.[^.]*/);
  $fn = $basename.".bed12";
  $bed12_fn = file($dest,$fn);
  open(BED12OUT, "> $bed12_fn");

  open(INBED, "< $bed_in") or croak $!;
  while(<INBED>){
    chomp;
    my ($chr, $start, $end, $info, $score, $strand) = split("\t");

    if ($info =~ /^splits\:(\d+)\:(\d+)\:(\d+):(\w):(\w)/){
      $reads = $1;
      $proper = $4;
      $passed = $5;
      if ($circ == 1){ # skip 'N' (normal), 'L' (left) and 'R' (right)
	next unless ($proper =~ /[C]$/);
      }
      else { # skip 'L', 'R' and 'C'
	next unless ($proper =~ /[N]$/);
      }
      next unless ($passed =~ /[PM]$/); # ignore 'F'
      next if($reads < $mcov);
    }
    else {
      croak "ERROR [$this_function] unsupported INFO field in input BED:\n$info\n";
    }
    $pos5 = $start;
    $pos3 = $end;

    @bedline = join("\t",$chr,eval($start-$window),
		    eval($end+$window),$info,$score,$strand,$start,$end,
		    "0","2","1,1","0,".eval($end-$start-1));
    print BED12OUT "@bedline\n";
  }
  close(INBED);
  close(BED12OUT);

  push (@result, $bed12_fn);
  return @result;
}

# intersect_sj($p_annot,$p_mapped,$dest,$prefix,$window,$mil)
#
# Intersect splice junctions determined by RNA-seq with annotated
# splice junctions. Determine novel and existing splice junctions.
#
# Writes BED6 files for existing and novel splice junctions to $dest
# and reurns an array with the absolute path to the two resulting BED
# files
sub intersect_sj{
  my ($p_annot,$p_mapped,$dest,$prefix,$window,$mil) = @_;
  my ($dha,$dhm);
  my $processed = 0;
  my @junctions = ();
  my @transcript_beds = ();
  my %asj = (); # annotated splice junctions hash
  my $bedtools = can_run('bedtools') or croak "bedtools not found";
  my $sortBed  = can_run('sortBed') or croak "sortBed not found";
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] $p_annot does not exist\n"
    unless (-d $p_annot);
  croak "ERROR [$this_function] $p_mapped does not exist\n"
    unless (-d $p_mapped);
  croak "ERROR [$this_function] $dest does not exist\n"
    unless (-d $dest);

  # get a list of all files in $p_annot
  opendir($dha, $p_annot) or croak "Cannot opendir $p_annot: $!";
  while(readdir $dha) { push @transcript_beds, $_; }
  closedir($dha);

  # get a list of all splice junctions seen in RNA-seq data
  opendir($dhm, $p_mapped) or croak "Cannot opendir $p_mapped: $!";
  @junctions = grep { /^(chr\d+)\_(\d+)-(\d+)\.mappedSS\.bed6/ } readdir($dhm);

  # process splice junctions seen in RNA-seq data
  foreach my $file (@junctions){
    $processed++;
    croak "Unexpected file name pattern\n" unless ($file =~ /(chr\d+)\_(\d+)-(\d+)/);
    my $sc = $1;
    my $s5 = $2;
    my $s3 = $3;
    my $pattern = $sc."_".$s5."-".$s3;

    my @annotated_beds = grep { /^(chr\d+)\_(\d+)-(\d+)/ && $2<=$s5 && $3>=$s3 && $sc eq $1} @transcript_beds;
    #print"\t intersecting against ".(eval $#annotated_beds+1)." transcripts: @annotated_beds \n";

    # intersect currently opened SJ against all transcripts in @annotated_beds
    foreach my $i (@annotated_beds){
      my $a = file($p_mapped,$file);
      my $b = file($p_annot,$i);
      my $intersect_cmd = "$bedtools intersect -a $a -b $b -c -nobuf";
      open(INTERSECT, $intersect_cmd."|");
      while(<INTERSECT>){
	chomp;
	if ($_ =~ /1$/) { $asj{$pattern} = 1;}
      }
      close(INTERSECT);
    }
    if ( $processed % 1000 == 0){
      print STDERR "processed $processed splice junctions\n";
    }
  }
  closedir($dhm);

  # go through the mapped splice junctions files once again and
  # separate novel from existing splice junctions
  if (length($prefix)>0){$prefix .= ".";}
  my $outname_exist   = file($dest, $prefix."exist.SS.bed");
  my $outname_exist_u = file($dest, $prefix."exist.SS.u.bed");
  my $outname_novel   = file($dest, $prefix."novel.SS.bed");
  my $outname_novel_u = file($dest, $prefix."novel.SS.u.bed");
  open (EXISTOUT, "> $outname_exist_u");
  open (NOVELOUT, "> $outname_novel_u");

  # write new ones to NOVELOUT; existing ones to EXISTOUT
  foreach my $file (@junctions){
    if ($file =~ m/^(chr\d+\_\d+-\d+)/){
       my $pattern = $1;
       my $fn = file($p_mapped,$file);
       open (SJ, "< $fn") or croak "Cannot open $fn $!";
       while(<SJ>){
	 chomp;
	 $_ = m/^(chr\w+)\s(\d+)\s(\d+)\s(splits:\d+:\d+:\d+:\w:\w)\s(\d+)\s([+-01])/;
	 my $chr = $1;
	 my $start = $2;
	 my $end = $3;
	 my $name = $4;
	 my $score = $5;
	 my $strand = $6;
	 my @bedline = join("\t",$chr,eval($start+$window),
			    eval($start+$window+1),$name,$score,$strand);
	 if (exists $asj{$pattern}){ # annotated splice junction
	   print EXISTOUT "@bedline\n";
	 }
	 else {  # novel splice junction
	   print NOVELOUT "@bedline\n";
	 }
       } # end while
       close(SJ);
     } # end if
    else{ carp "Error with parsing BED6 junction file names __FILE__ __LINE__\n";}
  }
  close(EXISTOUT);
  close(NOVELOUT);

  # sort the resulting bed files
  my $cmd = "$bedtools sort -i  $outname_exist_u > $outname_exist";
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $bedtools  unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  unlink($outname_exist_u);

  $cmd = "$bedtools sort -i  $outname_novel_u > $outname_novel";
  ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $bedtools  unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  unlink($outname_novel_u);

  printf STDERR "processed $processed splice junctions\n";
  my @result = ();
  push(@result, $outname_exist);
  push(@result, $outname_novel);
  return @result;
}

# ss_isCanonical ( $chr,$p5,$p3,$fastaO )

# Checks whether a given splice junction is canonical, ie. whether the
# first and last two nucleotides of the enclosed intron correspond to
# a certain nucleotide motif. $chr is the chromosome name, $p5 and $p3
# the 5' and 3' ends of the splice junction and $fastaobjR is a
# Bio::PrimarySeq::Fasta object holding the underlying reference
# genome.
#
# Th most common canonical splice junction motif is GT-AG (shown
# below). Other canonical motifs are GC->AG and AT->AC. 
#
#   ------------------>
# 5'===]GT..........AG[====3'
#
#   <-----------------
# 3'===]GA..........TG[====5'
#
sub ss_isCanonical{
  my ($chr,$p5,$p3,$fo) = @_;
  my ($seqL,$seqR,$pattern);
  my $ss_motif_length = 2;
  my $this_function = (caller(0))[3];
  my $c = -1;

  $seqL = $fo->stranded_subsequence($chr,$p5+1,$p5+$ss_motif_length,"+");
  $seqR = $fo->stranded_subsequence($chr,$p3-$ss_motif_length+1,$p3,"+");

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

Bio::ViennaNGS::SpliceJunc - Perl extension for alternative splicing
analysis

=head1 SYNOPSIS

  use Bio::ViennaNGS::SpliceJunc;
  use Bio::ViennaNGS::Fasta;

  # get a Bio::ViennaNGS::Fasta object
  my $fastaO = Bio::ViennaNGS::Fasta->new($fasta_in);

  # Extract annotated splice sites from BED12
  bed6_ss_from_bed12($bed12_in,$dest_annot,$window,$want_canonical,$fastaO);

  # Extract mapped splice junctions from RNA-seq data
  bed6_ss_from_rnaseq($bed_in,$dest_ss,$window,$mincov,$want_canonical,$fastaO);

  # Check for each splice junction seen in RNA-seq if it overlaps with
  # any annotated splice junction
  @res = intersect_sj($dest_annot,$dest_ss,$dest,$prefix,$window,$mil);

  # Convert splice junctions seen in RNA-seq data to BED12
  @res = bed6_ss_to_bed12($s_in,$outdir,$window,$mincov,$want_circular);

  # Check whether a splice junction is canonical
  $c = ss_isCanonical($chr,$pos5,$pos3,$fastaO);

=head1 DESCRIPTION

L<Bio::ViennaNGS::SpliceJunc> is a Perl module for alternative
splicing (AS) analysis. It provides routines for identification,
characterization and visualization of novel and existing (annotated)
splice junctions from RNA-seq data.

Identification of novel splice junctions is based on intersecting
potentially novel splice junctions from RNA-seq data with annotated
splice junctions.

=head2 SUBROUTINES

=over 3

=item bed6_ss_from_bed12($bed12,$dest,$window,$can,$fastaO)

Extracts splice junctions from a BED12 file (provided via argument
C<$bed12>), writes a BED6 file for each transcript to C<$dest>,
containing all its splice junctions. If C<$can> is 1, canonical splice
junctions are reported in the 'name' field of the output BED6 file.
Output splice junctions can be flanked by a window of +/- C<$window>
nt. C<$fastaO> is a L<Bio::ViennaNGS::Fasta> object. Each splice
junction is represented as two bed lines in the output BED6.

=item bed6_ss_from_rnaseq($bed_in,$dest,$window,$mcov,$can,$fastaO)

Extracts splice junctions from mapped RNA-seq data. The input BED6
file should contain coordinates of introns in the following syntax:

chr1    3913    3996    splits:97:97:97:N:P     0       +

The fourth column in this BED file (correponding to the 'name' field
according to the L<BED
specification|http://genome.ucsc.edu/FAQ/FAQformat.html#format1>)
should be a colon-separated string of six elements, where the first
element should be 'splits' and the second element is assumed to hold
the number of reads supporting this splice junction. The fifth element
indicates the splice junction type: A capital 'N' determines a normal
splice junction, whereas 'C' indicates circular and 'T' indicates
trans-splice junctions, respectively. Only normal splice junctions
('N') are considered, the rest is skipped. Elements 3, 4 and 6 are not
further processed.

We recommend using
F<segemehl|http://www.bioinf.uni-leipzig.de/Software/segemehl/> for
generating this type of BED6 files. This routine is, however, not
limited to F<segemehl> output. BED6 files containing splice junction
information from other short read mappers or third-party sources will
be processed if they are formatted as described above.

This routine writes a BED6 file for each splice junction provided in
the input to C<$dest>. Output splice junctions can be flanked by a
window of +/- C<$window> nt. Canonical splice junctions are reported
in the 'name' field of the output BED6 file if C<$can> is 1 and
C<$featO> is a L<Bio::ViennaNGS::Fasta> object. Each splice junction
is represented as two BED lines in the output BED6. Only splice
junctions that are supported by at least C<$mcov> reads are reported.

=item bed6_ss_to_bed12($bed_in,$dest,$window,$mcov,$circ)

Produce BED12 output for splice junctions found in RNA-seq data. Input
BED6 files (provided via C<$bed_in>) are supposed to conform to the
F<segemehl|http://www.bioinf.uni-leipzig.de/Software/segemehl/>
standard format for reporting splice junctions, which has the
following syntax:

chr1    3913    3996    splits:97:97:97:N:P     0       +

See L<bed6_ss_rom_rnaseq> for details.

C<$dest> is the output path. Output splice junctions can optionally be
flanked by a window of +/- C<$window> nt. Only splice junctions that
are supported by at least C<$mcov> reads are reported. If C<$circ> is
1, B<circular> splice junctions are reported (if present in the
input), else B<normal> splice junctions are processed.

=item intersect_sj($p_annot,$p_mapped,$dest,$prefix,$window,$mil)

Intersects all splice junctions identified in an RNA-seq experiment
with annotated splice junctions. Identifies and characterizes novel
and existing splice junctions. Each BED6 file in C<$p_mapped> is
intersected with those transcript splice junction BED6 files in
C<$p_annot>, whose genomic location spans the query splice
junction. This is to prevent the tool from intersecting each splice
site found in the mapped RNA-seq data with B<all> annotated
transcripts. C<$mil> specifies a maximum intron length.

The intersection operations are performed with F<bedtools intersect>
from the L<BEDtools|https://github.com/arq5x/bedtools2> suite). BED
sorting operations are performed with F<bedtools sort>.

Writes two BED6 files to C<$dest> (optionally prefixed by C<$prefix>), which
contain novel and existing splice junctions, respectively.

=item ss_isCanonical($chr,$p5,$p3,$fo)

Checks whether a given splice junction is canonical, ie. whether the
first and last two nucleotides of the enclosed intron correspond to a
certain nucleotide motif. C<$chr> is the chromosome name, C<$p5> and
C<$p3> the 5' and 3' ends of the splice junction and C<$fo> is a
L<Bio::ViennaNGS::Fasta> object holding the underlying reference
genome

This routine does not explicitly consider standedness in the sense
that splice junction motifs are evaluated in terms of the forward
strand of the underlying reference sequence. This is best explained by
an example: Consider the splice junction motif GU->G on the reverse
strand. In 5' to 3' direction of the forward strandm this junction
reads CT->AC. A splice junction is canonical if its motif corresponds
to one of the following cases:

5'===]GT|CT....AG|AC[====3' ie GT->AG or CT->AC
5'===]GC|CT....AG|GC[====3' ie GC->AG or CT->GC
5'===]AT|GT....AC|AT[====3' ie AT->AC or GT->AT

=back

=head1 DEPENDENCIES

This modules depends on the following Perl modules:

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::Fasta>

=item L<IPC::Cmd>

=item L<Path::Class>

=item L<Carp>

=back

L<Bio::ViennaNGS::SpliceJunc> uses third-party tools for computing
intersections of BED files: F<bedtools intersect> from the
L<BEDtools|http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html>
suite is used to compute overlaps and F<bedtools sort> is used to sort
BED output files. Make sure that those third-party utilities are
available on your system, and that hey can be found and executed by
the perl interpreter. We recommend installing the latest version of
L<BEDtools|https://github.com/arq5x/bedtools2> on your system.

=head1 SEE ALSO

L<Bio::ViennaNGS>

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.


=cut
