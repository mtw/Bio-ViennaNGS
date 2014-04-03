#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2014-04-03 12:55:16 mtw>
#
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

package ViennaNGS;

use Exporter;
use strict;
use warnings;
use Bio::Perl;
use Bio::DB::Sam;
use Data::Dumper;
use File::Basename qw(basename fileparse);
use File::Temp qw(tempfile);

our @ISA = qw(Exporter);
our $VERSION = '0.05';
our @EXPORT = qw(get_stranded_subsequence split_bam bam2bw bed2bw);

our @EXPORT_OK = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# get_stranded_subsequence ($obj,$start,$stop,$strand)
# retrieve RNA/DNA sequence from a Bio::PrimarySeqI /
# Bio::PrimarySeq::Fasta object
sub get_stranded_subsequence {
  my ($o,$start,$end,$strand) = @_;
  my $seq = $o->subseq($start => $end);
  if ($strand eq '-1' || $strand eq '-') {
    my $rc = revcom($seq);
    $seq = $rc->seq();
  }
  #print "id:$id\nstart:$start\nend:$end\n";
  return $seq;
}

# split_bam ( $bam,$reverse,$want_uniq,$want_bed,$dest_dir,$log )
# Splits BAM file $bam according to [+] and [-] strand
# Returns array with newly splitted BAM files
sub split_bam {
  my %data = ();
  my @NHval = ();
  my @processed_files = ();
  my $verbose = 0;
  my ($bamfile,$reverse,$want_uniq,$want_bed,$dest_dir,$log) = @_;
  my ($bam,$sam,$bn,$path,$ext,$header,$flag,$NH,$eff_strand,$tmp);
  my ($bam_pos,$bam_neg,$tmp_bam_pos,$tmp_bam_neg,$bamname_pos,$bamname_neg);
  my ($bed_pos,$bed_neg,$bedname_pos,$bedname_neg);
  my ($seq_id,$start,$stop,$strand,$target_names,$id,$score);
  my %count_entries = (
		       total     => 0,
		       uniq      => 0,
		       pos       => 0,
		       neg       => 0,
		       skip      => 0,
		       cur       => 0,
		       mult      => 0,
		       se_alis   => 0,
		       pe_alis   => 0,
		       flag      => 0,
		      );
  $data{count} = \%count_entries;
  $data{flag} = ();

  die("ERROR: $bamfile does not exist\n") unless (-e $bamfile);
  die("ERROR: $dest_dir does not exist\n") unless (-d $dest_dir);

  open(LOG, ">", $log) or die $!;

  (undef,$tmp_bam_pos) = tempfile('BAM_POS_XXXXXXX',UNLINK=>0);
  (undef,$tmp_bam_neg) = tempfile('BAM_NEG_XXXXXXX',UNLINK=>0);

  $bam = Bio::DB::Bam->open($bamfile, "r");
  $header = $bam->header;
  $target_names = $header->target_name;

  ($bn,$path,$ext) = fileparse($bamfile, qr /\..*/);
  unless ($dest_dir =~ /\/$/){$dest_dir .= "/";}
  $bamname_pos = $dest_dir.$bn.".pos".$ext;
  $bamname_neg = $dest_dir.$bn.".neg".$ext;
  $bam_pos = Bio::DB::Bam->open($tmp_bam_pos,'w')
    or die "Could not open bam_pos file for writing: $!";
  $bam_neg = Bio::DB::Bam->open($tmp_bam_neg,'w')
    or die "Could not open bam_neg file for writing: $!";

  if ($want_bed == 1){
     $bedname_pos = $dest_dir.$bn.".pos.bed";
     $bedname_neg = $dest_dir.$bn.".neg.bed";
     open($bed_pos, ">", $bedname_pos); open($bed_neg, ">", $bedname_neg);
  }

  $bam_pos->header_write($header);$bam_neg->header_write($header);
  if($reverse == 1) {  # switch +/- strand mapping
    $tmp = $bam_pos;$bam_pos = $bam_neg;$bam_neg = $tmp;
    $tmp = $bed_pos;$bed_pos = $bed_neg;$bed_neg = $tmp;
  }

  while (my $read= $bam->read1() ) {
    @NHval = ();
    $data{count}{total}++;
    if($verbose == 1){print STDERR $read->query->name."\t";}

    # check if NH (the SAM tag used to indicate multiple mappings) is set
    if ($read->has_tag("NH")) {
      @NHval = $read->get_tag_values("NH");
      $NH = $NHval[0];
      if ($NH == 1) {
	$data{count}{uniq}++;
	if ($verbose == 1) {print STDERR "NH:i:1\t";}
      }
      else {
	$data{count}{mult}++;
	if ($verbose == 1) {print STDERR "NH:i:".$NH."\t";}
	if ($want_uniq == 1) { # skip processing this read if it is a mutli-mapper
	  $data{count}{skip}++;
	  next;
	}
      }
      $data{count}{cur}++;
    }
    else{ warn "Read ".$read->query->name." does not have NH attribute\n";}

  #  print Dumper ($read->query);
    $strand = $read->strand;
    $seq_id = $target_names->[$read->tid];
    $start  = $read->start;
    $stop   = $read->end;
    $id     = "x"; # $read->qname; 
    $score  = 100;

    if ( $read->get_tag_values('PAIRED') ) { # paired-end
      if($verbose == 1) {print STDERR "pe\t";}
      $data{count}{pe_alis}++;
      if ( $read->get_tag_values('FIRST_MATE') ){ # 1st mate; take its strand as granted
	if($verbose == 1) {print STDERR "FIRST_MATE\t".$strand." ";}
	if ( $strand eq "1" ){
	  $bam_pos->write1($read);
	  if ($want_bed){printf $bed_pos "%s\t%d\t%d\t%s\t%d\t",$seq_id,eval($start-1),$stop,$id,$score;}
	  if ($reverse == 0){
	    $data{count}{pos}++; $eff_strand=$strand;
	    if ($want_bed){printf $bed_pos "%s\n", "+";}
	  }
	  else {
	    $data{count}{neg}++; $eff_strand=-1*$strand;
	    if ($want_bed){printf $bed_pos "%s\n", "-";}
	  }
	}
	elsif ($strand eq "-1") {
	  $bam_neg->write1($read);
	  if ($want_bed){printf $bed_neg "%s\t%d\t%d\t%s\t%d\t",$seq_id,eval($start-1),$stop,$id,$score;}
	  if ($reverse == 0){
	    $data{count}{neg}++; $eff_strand=$strand;
	    if ($want_bed){printf $bed_neg "%s\n", "-";}
	  }
	  else {
	    $data{count}{pos}++;$eff_strand=-1*$strand;
	    if ($want_bed){printf $bed_neg "%s\n", "+";}
	  }
	}
	else {die "Strand neither + nor - ...exiting!\n";}
      }
      else{ # 2nd mate; reverse strand since the fragment it belongs to is ALWAYS located
            # on the other strand
	if($verbose == 1) {print STDERR "SECOND_MATE\t".$strand." ";}
	if ( $strand eq "1" ) {
	  $bam_neg->write1($read);
	  if ($want_bed){printf $bed_neg "%s\t%d\t%d\t%s\t%d\t",$seq_id,eval($start-1),$stop,$id,$score;}
	  if ($reverse == 0){
	    $data{count}{neg}++;$eff_strand=$strand;
	    if ($want_bed){printf $bed_neg "%s\n", "-";}
	  }
	  else {
	    $data{count}{pos}++; $eff_strand=-1*$strand;
	    if ($want_bed){printf $bed_neg "%s\n", "+";}
	  }
	}
	elsif ( $strand eq "-1" ) {
	  $bam_pos->write1($read);
	  if ($want_bed){printf $bed_pos "%s\t%d\t%d\t%s\t%d\t",$seq_id,eval($start-1),$stop,$id,$score;}
	  if ($reverse == 0){
	    $data{count}{pos}++;$eff_strand=$strand;
	    if ($want_bed){printf $bed_pos "%s\n", "+";}
	  }
	  else {
	    $data{count}{neg}++;$eff_strand=-1*$strand;
	    if ($want_bed){printf $bed_pos "%s\n", "-";}
	  }
	}
	else {die "Strand neither + nor - ...exiting!\n";}
      }
    }
    else { # single-end
      if($verbose == 1) {print STDERR "se\t";}
      $data{count}{se_alis}++;
      if ( $read->strand eq "1" ){
	$bam_pos->write1($read);
	if ($want_bed){printf $bed_pos "%s\t%d\t%d\t%s\t%d\t",$seq_id,eval($start-1),$stop,$id,$score;}
	if ($reverse == 0){
	  $data{count}{pos}++;$eff_strand=$strand;
	  if ($want_bed){printf $bed_pos "%s\n", "+";}
	}
	else {
	  $data{count}{neg}++;$eff_strand=-1*$strand;
	  if ($want_bed){printf $bed_pos "%s\n", "-";}
	}
      }
      elsif ($read->strand eq "-1") {
	$bam_neg->write1($read);
	if ($want_bed){printf $bed_neg "%s\t%d\t%d\t%s\t%d\t",$seq_id,eval($start-1),$stop,$id,$score;}
	if ($reverse == 0){
	  $data{count}{neg}++;$eff_strand=$strand;
	  if ($want_bed){printf $bed_neg "%s\n", "-";}
	}
	else {
	  $data{count}{pos}++;$eff_strand=-1*$strand;
	  if ($want_bed){printf $bed_neg "%s\n", "+";}
	}
      }
      else {die "Strand neither + nor - ...exiting!\n";}
    }
    if($verbose == 1) {print STDERR "--> ".$eff_strand."\t";}

    # collect statistics of SAM flags
    $flag = $read->flag;
    unless (exists $data{flag}{$flag}){
      $data{flag}{$flag} = 0;
    }
    $data{flag}{$flag}++;
    if ($verbose == 1) {print STDERR "\n";}
  } # end while

  rename ($tmp_bam_pos, $bamname_pos);
  rename ($tmp_bam_neg, $bamname_neg);
  push (@processed_files, ($bamname_pos,$bamname_neg));
  if ($want_bed){
    push (@processed_files, ($bedname_pos,$bedname_neg))
  }

  # error checks
  unless ($data{count}{pe_alis} + $data{count}{se_alis} == $data{count}{cur}) {
    printf "ERROR:  paired-end + single-end alignments != total alignment count\n";
    print Dumper(\%data);
    die;
  }
  unless ($data{count}{pos} + $data{count}{neg} == $data{count}{cur}) {
    printf STDERR "%20d fragments on [+] strand\n",$data{count}{pos};
    printf STDERR "%20d fragments on [-] strand\n",$data{count}{neg};
    printf STDERR "%20d sum\n",eval($data{count}{pos}+$data{count}{neg});
    printf STDERR "%20d cur_count (should be)\n",$data{count}{cur};
    printf STDERR "ERROR: pos alignments + neg alignments != total alignments\n";
    print Dumper(\%data);
    die;
  }
  foreach (keys %{$data{flag}}){
    $data{count}{flag} += $data{flag}{$_};
  }
  unless ($data{count}{flag} == $data{count}{cur}){
    printf STDERR "%20d alignments considered\n",$data{count}{cur};
    printf STDERR "%20d alignments found in flag statistics\n",$data{count}{flag};
    printf STDERR "ERROR: #considered alignments != #alignments from flag stat\n";
    print Dumper(\%data);
    die;
  }

  # logging output
  printf LOG "# bam_split.pl log for file \'$bamfile\'\n";
  printf LOG "#-----------------------------------------------------------------\n";
  printf LOG "%20d total alignments (unique & multi-mapper)\n",$data{count}{total};
  printf LOG "%20d unique-mappers (%6.2f%% of total)\n",
    $data{count}{uniq},eval(100*$data{count}{uniq}/$data{count}{total}) ;
  printf LOG "%20d multi-mappers  (%6.2f%% of total)\n",
    $data{count}{mult},eval(100*$data{count}{mult}/$data{count}{total});
  printf LOG "%20d multi-mappers skipped\n", $data{count}{skip};
  printf LOG "%20d alignments considered\n", $data{count}{cur};
  printf LOG "%20d paired-end\n", $data{count}{pe_alis};
  printf LOG "%20s single-end\n", $data{count}{se_alis};
  printf LOG "%20d fragments on [+] strand  (%6.2f%% of considered)\n",
    $data{count}{pos},eval(100*$data{count}{pos}/$data{count}{cur});
  printf LOG "%20d fragments on [-] strand  (%6.2f%% of considered)\n",
    $data{count}{neg},eval(100*$data{count}{neg}/$data{count}{cur});
  printf LOG "#-----------------------------------------------------------------\n";
  printf LOG "Dumper output:\n". Dumper(\%data);
  close(LOG);
  return @processed_files;
}

# bam2bw ( $bam,$chromsizes )
# Generate BedGraph and BigWig coverage from BAM via two third-party tools:
# genomeCoverageBed from BEDtools
# bedGraphToBigWig from UCSC Genome Browser tools
sub bam2bw {
  my ($bamfile,$chromsizes) = @_;
  my $genomeCoverageBed = `which genomeCoverageBed`; chomp($genomeCoverageBed);
  my $bedGraphToBigWig = `which bedGraphToBigWig`; chomp($bedGraphToBigWig);
  my $outpath = "./";
  my $outfolder = $outpath."vis";
  my ($GCB_cmd,$BGTBW_cmd);

  unless (-e $bamfile) {
    die "ERROR: Cannot find $bamfile\n";
  }
  unless (-e $chromsizes) {
    die "ERROR: Cannot find $chromsizes ...BigWig cannot be generated\n";
  }
  mkdir $outfolder;

  $GCB_cmd = "$genomeCoverageBed -bg -ibam $bamfile -g $chromsizes > $outfolder/$bamfile.bg";
  $BGTBW_cmd = "$bedGraphToBigWig $outfolder/$bamfile.bg $chromsizes $outfolder/$bamfile.bw";

  print STDERR ">> $GCB_cmd\n>> $BGTBW_cmd\n";
  system($GCB_cmd);
  system($BGTBW_cmd);
}

# bed2bw ($bed,$chromsies,$strand,$dest_dir)
# Generate BedGraph and stranded BigWig coberage from BED via two third-party tools:
# genomeCoverageBed from BEDtools
# bedGraphToBigWig from UCSC Genome Browser tools
sub bed2bw {
  my ($bedfile,$chromsizes,$strand,$dest_dir) = @_;
  my ($bn,$path,$ext,$cmd);
  my $genomeCoverageBed = `which genomeCoverageBed`; chomp($genomeCoverageBed);
  my $bedGraphToBigWig = `which bedGraphToBigWig`; chomp($bedGraphToBigWig);
  my $awk = `which awk`; chomp($awk);

  print STDERR "##\$bedfile: $bedfile\n##\$chromsizes: $chromsizes\n##\$dest_dir: $dest_dir\n";
  unless (-e $bedfile) {
    die "ERROR: Cannot find $bedfile\n";
  }
  unless (-e $chromsizes) {
    die "ERROR: Cannot find $chromsizes ...BigWig cannot be generated\n";
  }
  unless (-d $dest_dir){
    die "ERROR: $dest_dir does not exist\n";
  }

  ($bn,$path,$ext) = fileparse($bedfile, qr /\..*/);

  if ($strand eq "+"){
    $cmd = "$genomeCoverageBed -bg -i $bedfile -g $chromsizes -strand $strand > $dest_dir/$bn.pos.bg";
    $cmd .= " && ";
    $cmd .= "$bedGraphToBigWig $dest_dir/$bn.pos.bg $chromsizes $dest_dir/$bn.pos.bw";
  }
  else{
    $cmd = "$genomeCoverageBed -bg -i $bedfile -g $chromsizes -strand $strand > $dest_dir/$bn.neg.bg.1";
    $cmd .= " && cat $dest_dir/$bn.neg.bg.1 | $awk \'{ \$4 = - \$4 ; print \$0 }\' > $dest_dir/$bn.neg.bg";
    $cmd .= " && $bedGraphToBigWig $dest_dir/$bn.neg.bg $chromsizes  $dest_dir/$bn.neg.bw";
    unlink("$dest_dir/$bn.neg.bg.1");
  }
  print STDERR ">>$cmd\n";
  system($cmd);
}

1;
__END__


=head1 NAME

ViennaNGS - Perl extension for analysis of Next-Generation Sequencing
(NGS) data.

=head1 SYNOPSIS

  use ViennaNGS;

=head1 DESCRIPTION

ViennaNGS is a collection of subroutines often used for NGS data analysis.

=head1 EXPORT

Routines: get_stranded_subsequence($obj,$start,$stop,$path)
          split_bam($bam,$reverse,$want_uniq,$dest_dir,$log)
          bam2bw($bam,$chromsizes)
          bed2bw($bed,$chromsizes,$strand,$dest_dir)

Variables: none

=head2 get_stranded_subsequence($object,$start,$stop,$strand)

Returns the actual DNA/RNA sequence from $start to $stop. $object is a
Bio::PrimarySeq::Fasta object, which obeys the Bio::PrimarySeqI
conventions. To recover the entire raw DNA or protein sequence,
call $object->seq(). $strand is 1 or -1.

=head2 split_bam($bam,$reverse,$want_uniq,$dest_dir,$log)

Splits BAM file $bam according to [+] and [-] strand. $reverse,
$want_uniq and $want_bed are switches with values of 0 or 1,
triggering forced reversion of strand mapping (due to RNA-seq protocol
constraints), filtering of unique mappers (identified via NH:i:1 SAM
argument), and forced output of a BED file corresponding to
strand-specific mapping, respectively. $log holds name and path of the
log file.

Strand-splitting is done in a way that in paired-end alignments, FIRST
and SECOND mates (reads) are treated as _one_ fragment, ie FIRST_MATE
reads determine the strand, while SECOND_MATE reads are assigned the
opposite strand per definitionem. (This also holds if the reads are
not mapped in proper pairs and even if ther is no mapping partner at
all)

Sometimes the library preparation protocol causes inversion of the
read assignment (with respect to the underlying annotation). In those
cases, the natural mapping of the reads can be obtained by the
$reverse flag.

NOTE: Filtering of unique mappers is only safe for single-end
experiments; In paired-end experiments, read and mate are treated
separately, thus allowing for scenarios where eg. one read is a
multi-mapper, whereas its associate mate is a unique mapper, resulting
in an ambiguous alignment of the entire fragment.

=head2 bam2bw($bam,$chromsizes)

Creates BedGraph and BigWig coverage files from BAM. These can easily
be visualized as TrackHubs within the UCSC Genome Browser (see
http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html). Internally,
the conversion is accomplished by two third-party applications:
genomeCoverageBed (from BEDtools, see
http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html)
and bedGraphToBigWig (from the UCSC Genome Browser, see
http://hgdownload.cse.ucsc.edu/admin/exe/).

=head2 bed2bw($bed,$chromsizes,$strand,$dest_dir)

Creates BedGraph and stranded BigWig coverage profiles from BED
files. $chromsizes is the chromosome.sizes files, $strand is either
"+" or "-", and $dest_dir contains the output path for results.

Stranded BigWigs can easily be visualized via TrackHubs in the
UCSC Genome Browser. Internally, the conversion is accomplished by two
third-party applications: genomeCoverageBed (from BEDtools, see
http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html)
and bedGraphToBigWig (from the UCSC Genome Browser, see
http://hgdownload.cse.ucsc.edu/admin/exe/).

=head1 SEE ALSO

perldoc ViennaNGS::AnnoC

=head1 AUTHOR

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Michael Thomas Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.


=cut
