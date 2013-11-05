#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2013-11-05 12:02:12 mtw>
#
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael Thomas Wolfinger <michael@wolfinger.eu>
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
our $VERSION = '0.02';
our @EXPORT = qw(get_stranded_subsequence split_bam);

our @EXPORT_OK = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# get_stranded_subsequence ($obj,$start,$stop,$path)
# retrieve RNA/DNA sequence from a Bio::PrimarySeqI /
# Bio::PrimarySeq::Fasta object
sub get_stranded_subsequence {
  my ($o,$start,$end,$strand) = @_;
  my $seq = $o->subseq($start => $end);
  if ($strand == -1) {
    my $rc = revcom($seq);
    $seq = $rc->seq();
  }
  #print "id:$id\nstart:$start\nend:$end\n";
  return $seq;
}

# split_bam ( $bam,$reverse,$want_uniq,$log )
# Splits BAM file $bam according to [+] and [-] strand
# Returns array with newly splitted BAM files
sub split_bam {
  my %data = ();
  my @processed_bam = ();
  my ($bamfile,$reverse,$want_uniq,$log) = @_;
  my ($bam,$sam,$bn,$path,$ext,$header,$flag);
  my ($bam_pos,$bam_neg,$tmpname_pos,$tmpname_neg,$bamname_pos,$bamname_neg);
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

  open(LOG, ">", $log) or die $!;

  (undef,$tmpname_pos) = tempfile('BAM_POS_XXXXXXX',UNLINK=>0);
  (undef,$tmpname_neg) = tempfile('BAM_NEG_XXXXXXX',UNLINK=>0);

  $bam = Bio::DB::Bam->open($bamfile, "r");
  $header = $bam->header;

  ($bn,$path,$ext) = fileparse($bamfile, qr /\..*/);
  $bamname_pos = $bn.".pos".$ext; $bamname_neg = $bn.".neg".$ext;
  $bam_pos = Bio::DB::Bam->open($tmpname_pos,'w')
    or die "Could not open bam_pos file for writing: $!";
  $bam_neg = Bio::DB::Bam->open($tmpname_neg,'w')
    or die "Could not open bam_neg file for writing: $!";

  $bam_pos->header_write($header);$bam_neg->header_write($header);

  while (my $read= $bam->read1() ) {
    my @NHval = ();
    my $NH;
    $data{count}{total}++;

    # check if NH (the SAM tag used to indicate multiple mappings) is set
    if ($read->has_tag("NH")) {
      @NHval = $read->get_tag_values("NH");
      $NH = $NHval[0];
      if ($NH == 1) {$data{count}{uniq}++;}
      else {
	$data{count}{mult}++;
	if ($want_uniq == 1) { # skip processing this read if it is a mutli-mapper
	  $data{count}{skip}++;
	  next;
	}
      }
      $data{count}{cur}++;
      # print $aln->query->name $NH \n";
    }
    else{
      warn "Read ".$read->query->name." does not have NH attribute\n";
    }

    if($reverse == 1) {  # switch +/- strand mapping
      my $tmp = $bam_pos;$bam_pos = $bam_neg;$bam_neg = $tmp;
    }

    if ( $read->get_tag_values('PAIRED') ) { # paired-end
      $data{count}{pe_alis}++;
      if ( $read->get_tag_values('FIRST_MATE') ){ # 1st mate; take its strand as granted
	if ( $read->strand eq "1" ){
	  $bam_pos->write1($read);
	  ($reverse == 0) ? $data{count}{pos}++ : $data{count}{neg}++;
	}
	elsif ($read->strand eq "-1") {
	  $bam_neg->write1($read);
	  ($reverse == 0) ? $data{count}{neg}++ : $data{count}{pos}++;
	}
	else {die "Strand neither + nor - ...exiting!\n";}
      }
      else{ # 2nd mate; reverse strand since the fragment it belongs to is ALWAYS located
            # on the other strand
	if ( $read->strand eq "1" ) {
	  $bam_neg->write1($read);
	  ($reverse == 0) ? $data{count}{neg}++ : $data{count}{pos}++;
	}
	elsif ( $read->strand eq "-1" ) {
	  $bam_pos->write1($read);
	  ($reverse == 0) ? $data{count}{pos}++ : $data{count}{neg}++
	}
	else {die "Strand neither + nor - ...exiting!\n";}
      }
    }
    else { # single-end
      $data{count}{se_alis}++;
      if ( $read->strand eq "1" ){
	$bam_pos->write1($read);
	($reverse == 0) ? $data{count}{pos}++ : $data{count}{neg}++;
      }
      elsif ($read->strand eq "-1") {
	$bam_neg->write1($read);
	($reverse == 0) ? $data{count}{neg}++ : $data{count}{pos}++;
      }
      else {die "Strand neither + nor - ...exiting!\n";}
    }

    # collect statistics of SAM flags
   $flag = $read->flag;
   unless (exists $data{flag}{$flag}){
      $data{flag}{$flag} = 0;
    }
    $data{flag}{$flag}++;

  } # end while

  rename ($tmpname_pos, $bamname_pos);
  rename ($tmpname_neg, $bamname_neg);
  push (@processed_bam, ($bamname_pos,$bamname_neg));

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
  return @processed_bam;
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
          split_bam($bam,$reverse,$want_uniq,$log)

Variables: none

=head2 get_stranded_subsequence($object,$start,$stop,$strand)

Returns the actual DNA/RNA sequence from $start to $stop. $object is a
Bio::PrimarySeq::Fasta object, which obeys the Bio::PrimarySeqI
conventions. To recover the entire raw DNA or protein sequence,
call $object->seq(). $strand is 1 or -1.

=head2 split_bam($bam,$reverse,$want_uniq,$log)

Splits BAM file $bam according to [+] and [-] strand. $reverse and
$want_uniq are switches with values of 0 or 1, triggering forced
reversion of strand mapping (due to RNA-seq protocol constraints) and
filtering os unique mappers (identified via NH:i:1 SAM argument),
respectively. NOTE: Filtering of unique mappers is only safe for
single-end experiments; In paired-end experiments, read and mate are
treated separately, thus allowing for scenarios where eg. one read is
a multi-mapper, whereas its associate mate is a unique mapper,
resulting in an ambiguous alignment of the entire fragment. $log hold
the name (and path) of the log file.

=head1 SEE ALSO

perldoc ViennaNGS::AnnoC

=head1 AUTHOR

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Michael Thomas Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.


=cut
