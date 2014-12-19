# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:31:03 mtw>

package Bio::ViennaNGS::Bam;

use 5.12.0;
use Exporter;
use version; our $VERSION = qv('0.12_07');
use strict;
use warnings;
use Bio::Perl 1.00690001;
use Bio::DB::Sam 1.39;
use Data::Dumper;
use File::Basename qw(fileparse);
use File::Temp qw(tempfile);
use Path::Class;
use Carp;

our @ISA = qw(Exporter);
our @EXPORT = ();

our @EXPORT_OK = qw ( split_bam uniquify_bam );

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

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
  my $this_function = (caller(0))[3];
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

  croak "ERROR [$this_function] $bamfile does not exist\n"
    unless (-e $bamfile);
  croak "ERROR [$this_function] $dest_dir does not exist\n"
    unless (-d $dest_dir);

  open(LOG, ">", $log) or croak $!;

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
    or croak "ERROR [$this_function] Could not open bam_pos file for writing: $!";
  $bam_neg = Bio::DB::Bam->open($tmp_bam_neg,'w')
    or croak "ERROR [$this_function] Could not open bam_neg file for writing: $!";

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
    else{ carp "WARN [$this_function] Read ".$read->query->name." does not have NH attribute\n";}

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
	else {croak "Strand neither + nor - ...exiting!\n";}
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
	else {croak "Strand neither + nor - ...exiting!\n";}
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
      else {croak "Strand neither + nor - ...exiting!\n";}
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
  push (@processed_files, ($data{count}{pos},$data{count}{neg}));
  if ($want_bed){
    push (@processed_files, ($bedname_pos,$bedname_neg))
  }

  # error checks
  unless ($data{count}{pe_alis} + $data{count}{se_alis} == $data{count}{cur}) {
    printf "ERROR:  paired-end + single-end alignments != total alignment count\n";
    print Dumper(\%data);
    croak $!;
  }
  unless ($data{count}{pos} + $data{count}{neg} == $data{count}{cur}) {
    printf STDERR "%20d fragments on [+] strand\n",$data{count}{pos};
    printf STDERR "%20d fragments on [-] strand\n",$data{count}{neg};
    printf STDERR "%20d sum\n",eval($data{count}{pos}+$data{count}{neg});
    printf STDERR "%20d cur_count (should be)\n",$data{count}{cur};
    printf STDERR "ERROR: pos alignments + neg alignments != total alignments\n";
    print Dumper(\%data);
    croak $!;
  }
  foreach (keys %{$data{flag}}){
    $data{count}{flag} += $data{flag}{$_};
  }
  unless ($data{count}{flag} == $data{count}{cur}){
    printf STDERR "%20d alignments considered\n",$data{count}{cur};
    printf STDERR "%20d alignments found in flag statistics\n",$data{count}{flag};
    printf STDERR "ERROR: #considered alignments != #alignments from flag stat\n";
    print Dumper(\%data);
    croak $!;
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

sub uniquify_bam {
  my ($bamfile,$dest,$log) = @_;
  my ($bam, $bn,$path,$ext,$read,$header);
  my ($tmp_uniq,$tmp_mult,$fn_uniq,$fn_mult,$bam_uniq,$bam_mult);
  my ($count_all,$count_uniq,$count_mult) = (0)x3;
  my @processed_files = ();
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] Cannot find $bamfile\n"
    unless (-e $bamfile);
  croak "ERROR [$this_function] $dest does not exist\n"
    unless (-d $dest);

  ($bn,$path,$ext) = fileparse($bamfile, qr /\..*/);

  (undef,$tmp_uniq) = tempfile('BAM_UNIQ_XXXXXXX',UNLINK=>0);
  (undef,$tmp_mult) = tempfile('BAM_MULT_XXXXXXX',UNLINK=>0);

  $bam     = Bio::DB::Bam->open($bamfile, "r");
  $fn_uniq = file($dest,$bn.".uniq".$ext);
  $fn_mult = file($dest,$bn.".mult".$ext);
  $header  = $bam->header; # TODO: modify header, leave traces ...

  $bam_uniq = Bio::DB::Bam->open($tmp_uniq,'w')
    or croak "ERROR [$this_function] Cannot open temp file for writing: $!";
  $bam_mult = Bio::DB::Bam->open($tmp_mult,'w')
    or croak "ERROR [$this_function] Cannot open temp file for writing: $!";
  $bam_uniq->header_write($header);
  $bam_mult->header_write($header);

  while ($read = $bam->read1() ) {
    $count_all++;
    if ($read->aux_get("NH") == 1){ # uniquely mapped reads
      $bam_uniq->write1($read);
      $count_uniq++;
    }
    else { # multiply mapped reads
      $bam_mult->write1($read);
      $count_mult++;
    }
  }

  croak "ERROR [$this_function] Read counts don't match\n"
    unless ($count_uniq + $count_mult == $count_all);

  rename ($tmp_uniq, $fn_uniq);
  rename ($tmp_mult, $fn_mult);
  push (@processed_files, ($fn_uniq,$fn_mult));

  if (defined $log){
    my $lf = file($dest,$log);
    open(LOG, ">>", $lf) or croak $!;
    printf LOG "%15d reads total\n%15d unique reads\n%15d multiple reads\n",
      $count_all,$count_uniq,$count_mult;
    close(LOG);
  }
}


1;
__END__


=head1 NAME

Bio::ViennaNGS::Bam - High-level access to BAM files

=head1 SYNOPSIS

  use Bio::ViennaNGS::Bam;

  # split a single-end  or paired-end BAM file by strands
  @result = split_bam($bam_in,$rev,$want_uniq,$want_bed,$destdir,$logfile);

  # extract unique and multi mappers from a BAM file
  @result = uniquify_bam($bam_in,$outdir,$logfile);

=head1 DESCRIPTION

Bio::ViennaNGS::BAM provides high-level access to BAM file. Building
on L<Bio::DB::Sam>, it provides code to extract specific portions from
BAM files. It comes with routines for splitting BAM files by strand
(which is often required for visualization of NGS data) and extracting
uniquely and multiply aligned reads from BAM files.

=head2 ROUTINES

=over

=item split_bam($bam,$reverse,$want_uniq,$want_bed,$dest_dir,$log)

Splits BAM file $bam according to [+] and [-] strand. C<$reverse>,
C<$want_uniq> and C<$want_bed> are switches with values of 0 or 1,
triggering forced reversion of strand mapping (due to RNA-seq protocol
constraints), filtering of unique mappers (identified via NH:i:1 SAM
argument), and forced output of a BED file corresponding to
strand-specific mapping, respectively. C<$log> holds name and path of
the log file.

Strand-splitting is done in a way that in paired-end alignments, FIRST
and SECOND mates (reads) are treated as _one_ fragment, ie FIRST_MATE
reads determine the strand, while SECOND_MATE reads are assigned the
opposite strand I<per definitionem>. This also holds if the reads are
not mapped in proper pairs and even if there is no mapping partner at
all.

Sometimes the library preparation protocol causes inversion of the
read assignment (with respect to the underlying annotation). In those
cases, the natural mapping of the reads can be obtained by the
C<$reverse> flag.

This routine returns an array whose fist two elements are the file
names of the newly generate BAM files with reads mapped to the
positive, and negative strand, respectively. Elements three and four
are the number of fragments mapped to the positive and negative
strand. If the C<$want_bed> option was given elements fiveand six are
the file names of the output BED files for positive and negative
strand, respectively.

NOTE: Filtering of unique mappers is only safe for single-end
experiments; In paired-end experiments, read and mate are treated
separately, thus allowing for scenarios where eg. one read is a
multi-mapper, whereas its associate mate is a unique mapper, resulting
in an ambiguous alignment of the entire fragment.

=item uniquify_bam($bam,$dest,$log)

Extract I<uniquely> and I<multiply> aligned reads from BAM file
C<$bam> by means of the I<NH:i:> SAM attribute. New BAM files for
unique and multi mappers are created in the output folder C<$dest>,
which are named B<basename.uniq.bam> and B<basename.mult.bam>,
respectively. If defined, a logfile named C<$log> is created in the
output folder.

This routine returns an array holding file names of the newly created
BAM files for unique and multi mappers, respectively.

NOTE: Not all short read mappers use the I<NH:i:> SAM attribute to
decorate unique and multi mappers. As such, this routine will not work
unless your BAM file has these attributes.

=head1 DEPENDENCIES

=over 6

=item  L<Bio::Perl> >= 1.00690001

=item  L<BIO::DB::Sam> >= 1.39

=item  L<File::Basename>

=item  L<File::Temp>

=item  L<Path::Class>

=item  L<Carp>

=back

=head1 AUTHORS

=over

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
