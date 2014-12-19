# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:34:19 mtw>

package Bio::ViennaNGS::Util;

use 5.12.0;
use Exporter;
use version; our $VERSION = qv('0.12_07');
use strict;
use warnings;
use Bio::Perl 1.00690001;
use Bio::DB::Sam 1.39;
use Data::Dumper;
use File::Basename qw(basename fileparse);
use File::Temp qw(tempfile);
use IPC::Cmd qw(can_run run);
use Path::Class;
use Math::Round;
use Carp;
use Bio::ViennaNGS::FeatureChain;

our @ISA = qw(Exporter);
our @EXPORT = ();

our @EXPORT_OK = qw ( bed_or_bam2bw sortbed bed2bigBed computeTPM
		      featCount_data parse_multicov write_multicov
		      unique_array kmer_enrichment extend_chain
		      parse_bed6 fetch_chrom_sizes);


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

our @featCount = ();
my %unique = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub featCount_data {
  return \@featCount;
}

sub bed_or_bam2bw {
  my ($type,$infile,$chromsizes,$strand,$dest,$want_norm,$size,$scale,$log) = @_;
  my ($fn_bg_tmp,$fn_bg,$fn_bw);
  my ($bn,$path,$ext,$cmd);
  my @processed_files = ();
  my $factor = 1.;
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] \$type is '$type', however it is expected to be either 'bam' or 'bed'\n"
    unless ($type eq "bam") || ($type eq "bed");

  my $genomeCoverageBed = can_run('genomeCoverageBed') or
    croak "ERROR [$this_function] genomeCoverageBed utility not found";
  my $bedGraphToBigWig = can_run('bedGraphToBigWig') or
    croak "ERROR [$this_function] bedGraphToBigWig utility not found";
  my $awk = can_run('awk') or
    croak "ERROR [$this_function] awk utility not found";

  if(defined $log){
    open(LOG, ">>", $log) or croak $!;
    print LOG "LOG [$this_function] \$infile: $infile\n";
    print LOG "LOG [$this_function] \$dest: $dest\n";
    print LOG "LOG [$this_function] \$chromsizes: $chromsizes\n";
  }

  croak "ERROR [$this_function] Cannot find $infile\n"
    unless (-e $infile);
  croak "ERROR [$this_function] $dest does not exist\n"
    unless (-d $dest);
  croak "ERROR [$this_function] Cannot find $chromsizes\n"
      unless (-e $chromsizes);

  if ($want_norm == 1){
    $factor = $scale/$size;
    print LOG "LOG [$this_function] normalization: $factor = ($scale/$size)\n"
      if(defined $log);
  }

  ($bn,$path,$ext) = fileparse($infile, qr /\..*/);
  $fn_bg_tmp  = file($dest,$bn.".tmp.bg");
  $fn_bg      = file($dest,$bn.".bg");
  if($strand eq "+"){
    $fn_bw  = file($dest,$bn.".pos.bw");
  }
  else {
    $fn_bw  = file($dest,$bn.".neg.bw");
  }

  $cmd = "$genomeCoverageBed -bg -scale $factor -split ";
  if ($type eq "bed"){ $cmd .= "-i $infile -g $chromsizes"; } # chrom.sizes only required for processing BED
  else { $cmd .= "-ibam $infile "; }
  $cmd .= " > $fn_bg_tmp";

  if($strand eq "-"){
    $cmd .= " && cat $fn_bg_tmp | $awk \'{ \$4 = - \$4 ; print \$0 }\' > $fn_bg";
  }
  else{
    $fn_bg = $fn_bg_tmp;
  }
  $cmd .= " && $bedGraphToBigWig $fn_bg $chromsizes $fn_bw";

  if (defined $log){ print LOG "LOG [$this_function] $cmd\n";}

  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );

  if( !$success ) {
    print STDERR "ERROR [$this_function] External command call unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }
  if (defined $log){ close(LOG); }

  unlink ($fn_bg_tmp);
  unlink ($fn_bg);
  return $fn_bw;
}

sub bed2bigBed {
  my ($infile,$chromsizes,$dest,$log) = @_;
  my ($bn,$path,$ext,$cmd,$outfile);
  my $this_function = (caller(0))[3];
  my $bed2bigBed = can_run('bedToBigBed') or
    croak "ERROR [$this_function] bedToBigBed utility not found";

  if (defined $log){
    open(LOG, ">>", $log) or croak $!;
    print LOG "LOG [$this_function] \$infile: $infile -- \$chromsizes: $chromsizes --\$dest: $dest\n";
  }

  croak "ERROR [$this_function] Cannot find $infile"
    unless (-e $infile);
  croak "ERROR [$this_function] Cannot find $chromsizes"
    unless (-e $chromsizes);
  croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);

  # .bed6 .bed12 extensions are replaced by .bb
  ($bn,$path,$ext) = fileparse($infile, qr /\.bed[126]?/);
  $outfile = file($dest, "$bn.bb");

  $cmd = "$bed2bigBed $infile -extraIndex=name -tab $chromsizes $outfile";
  if (defined $log){ print LOG "LOG [$this_function] $cmd\n"; }
  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );

  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $bed2bigBed unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }

  if (defined $log){ close(LOG); }

  return $outfile;
}

sub sortbed {
  my ($infile,$dest,$outfile,$rm_orig,$log) = @_;
  my ($cmd,$out);
  my $this_function = (caller(0))[3];
  my $bedtools = can_run('bedtools') or
    croak "ERROR [$this_function bedtools utility not found";

  croak "ERROR [$this_function] Cannot find $infile"
    unless (-e $infile);
  croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);
  if (defined $log){open(LOG, ">>", $log) or croak $!;}

  $out = file($dest,$outfile);
  $cmd = "$bedtools sort -i $infile > $out";
  if (defined $log){ print LOG "LOG [$this_function] $cmd\n"; }
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] Call to $bedtools unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }

  if($rm_orig){
    unlink($infile) or
      carp "WARN [$this_function] Could not unlink $infile";
    if (defined $log){
      print LOG "[$this_function] removed $infile $!\n";
    }
  }

  if (defined $log){ close(LOG); }
}

sub computeTPM {
  my ($featCount_sample,$rl) = @_;
  my ($TPM,$T,$totalTPM) = (0)x3;
  my ($i,$features,$meanTPM);

  $features = keys %$featCount_sample; # of of features in hash
  print Dumper(\%$featCount_sample);print ">>$rl<<\n";

  # iterate through $featCount_sample twice:
  # 1. for computing T (denominator in TPM formula)
  foreach $i (keys %$featCount_sample){
    $T += ($$featCount_sample{$i}{count} * $rl)/($$featCount_sample{$i}{len});
  }
  # 2. for computng actual TPM values
  foreach $i (keys %$featCount_sample){
    $TPM = 1000000 * $$featCount_sample{$i}{count} * $rl/($$featCount_sample{$i}{len} * $T);
    $$featCount_sample{$i}{TPM} = $TPM;
    $totalTPM += $TPM;
  }
  $meanTPM = $totalTPM/$features;
  # print "totalTPM=$totalTPM | meanTPM=$meanTPM\n";
  return $meanTPM;
}

sub parse_multicov {
  my ($file) = @_;
  my @mcData = ();
  my ($mcSamples,$i);
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] multicov file $file not available\n"
    unless (-e $file);
  open (MULTICOV_IN, "< $file") or croak $!;

  while (<MULTICOV_IN>){
    chomp;
    @mcData = split(/\t/); # 0:chr|1:start|2:end|3:name|4:score|5:strand
    $mcSamples = (scalar @mcData)-6; # multicov extends BED6
    #print "$_\n";
    for ($i=0;$i<$mcSamples;$i++){
      $featCount[$i]{$mcData[3]} = {
				    chr    => $mcData[0],
				    start  => $mcData[1],
				    end    => $mcData[2],
				    name   => $mcData[3],
				    score  => $mcData[4],
				    strand => $mcData[5],
				    len    => eval($mcData[2]-$mcData[1]),
				    count  => $mcData[eval(6+$i)],
				   }
    }
    #print Dumper(@mcData);
  }
  close(MULTICOV_IN);
  return $mcSamples;
}

sub write_multicov {
  my ($item,$dest,$base_name) = @_;
  my ($outfile,$mcSamples,$nrFeatures,$feat,$i);
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function]: $dest does not exist\n"
    unless (-d $dest);
  $outfile = file($dest,$base_name.".".$item.".multicov.csv");
  open (MULTICOV_OUT, "> $outfile") or croak $!;

  $mcSamples = scalar @featCount; # of samples in %{$featCount}
  $nrFeatures = scalar keys %{$featCount[1]}; # of keys in %{$featCount}[1]
  #print "=====> write_multicov: writing multicov file $outfile with $nrFeatures lines and $mcSamples conditions\n";

  # check whether each column in %$featCount has the same number of entries
  for($i=0;$i<$mcSamples;$i++){
    my $fc = scalar keys %{$featCount[$i]}; # of keys in %{$featCount}
    #print "condition $i => $fc keys\n";
    unless($nrFeatures == $fc){
      croak "ERROR [$this_function]: unequal element count in \%\$featCount\nExpected $nrFeatures have $fc in condition $i\n";
    }
  }

  foreach $feat (keys  %{$featCount[1]}){
    my @mcLine = ();
    # process standard BED6 fields first
    push @mcLine, (${$featCount[1]}{$feat}->{chr},
		   ${$featCount[1]}{$feat}->{start},
		   ${$featCount[1]}{$feat}->{end},
		   ${$featCount[1]}{$feat}->{name},
		   ${$featCount[1]}{$feat}->{score},
		   ${$featCount[1]}{$feat}->{strand});
    # process multicov values for all samples

    for($i=0;$i<$mcSamples;$i++){
     # print "------------>  ";  print "processing $i th condition ";  print "<-----------\n";
      unless (defined ${$featCount[$i]}{$feat}){
	croak "Could not find item $feat in mcSample $i\n";
      }
      push @mcLine, ${$featCount[$i]}{$feat}->{$item};

    }
    #print Dumper(\@mcLine);
    print MULTICOV_OUT join("\t",@mcLine)."\n";
  }
  close(MULTICOV_OUT);
}

sub unique_array{

    my $arrayref = shift;
    my @array = @{$arrayref};

    foreach my $item (@array)
    {
	$unique{$item} ++;
    }
    my @arrayuid = sort {$a cmp $b} keys %unique;

    return(\@arrayuid);
}

sub kmer_enrichment{

    my @seqs =  @{$_[0]};
    my $klen     = $_[1]; 
#    my @seq = split( //, $read_tmp );
    my $kstring ='';
#return variables
    my %km;
    foreach my $sequences (@seqs){
#      print STDERR $sequences,"\n";
      my @seq = split( //, $sequences );
      for ( my $seq_pos = 0; $seq_pos <= $#seq-$klen ; $seq_pos++ ) {
	for (my $i=$seq_pos;$i<=$seq_pos+($klen-1);$i++){
	  $kstring .= $seq[$i]; 
	}
	$km{$kstring}++;
	$kstring = "";
      }
    }
    return( \%km );
}

sub extend_chain{
  my %sizes = %{$_[0]};
  my $chain = $_[1];
  my $l	    = $_[2];
  my $r	    = $_[3];
  my $u     = $_[4];
  my $d     = $_[5];

  ##return a new chain with extended coordinates
  my $extendchain = $chain -> clone();
  ## got through all features in original chain, calculate new start and end and safe in extendchain
  my @featarray = @{$extendchain->chain};
  foreach my $feature (@featarray){
    my $chrom  = $feature->chromosome;
    my $start  = $feature->start;
    my $end    = $feature->end;
    my $strand = $feature->strand;
    my $right  = 0;
    my $left   = 0;
    my $width  = nearest(1,($end-$start)/2);
    $width = 0 if ($d > 0 || $u > 0);
    if ($strand eq "+"){
      if ($d > 0){
	$start = $end;
	$r = $d;
      }
      if ($u > 0){
	$end = $start;
	$l = $u;
      }
      $right=$r;
      $left=$l;
    }
    elsif ($strand eq "-"){
      if ($u > 0){
	$start = $end;
	$l = $u;
      }
      if ($d > 0){
	$end = $start;
	$r = $d;
      }
      $right=$l;
      $left=$r;
    }
    if (($right-$width) <= 0){
      $right = 0;
    }
    else{
      $right-=$width;
    }
    if (($left-$width) <= 0 ){
      $left = 0;
    }
    else{
      $left-=$width;
    }
    if ( $start-$left >= 1 ){
      if ($end+$right >= $sizes{"chr".$chrom}){
	$end = $sizes{"chr".$chrom};
      }
      else{
	$end += $right;
      }
      $start -= $left;
    }
    elsif ( $start-$left <= 0 ){
      $start = 0;
      if ($end+$right >= $sizes{"chr".$chrom}){
	$end = $sizes{"chr".$chrom};
      }
      else{
	$end = $end+$right;
      }
    }
    else{
      die "Something wrong here!\n";
    }
    $feature->start($start);
    $feature->end($end);
  }
  $extendchain->type('extended');
  return($extendchain);
}

sub parse_bed6{
  my $bedfile = shift;
  open (my $Bed, "<:gzip(autopop)",$bedfile) or die "$!";
  my @featurelist; ## This will become a FeatureChain
  while(<$Bed>){
    ### This should be done by FeatureIO
    chomp (my $raw = $_);
    push my @line , split (/\t/,$raw);
    push @line, "\." if ( !$line[5] ); 

    (my $chromosome  = $line[0])=~ s/chr//g;
    my $start	     = $line[1];
    my $end	     = $line[2];
    my $name	     = $line[3];
    my $score	     = $line[4];
    my $strand	     = $line[5];
    my $extension = '';

    if ($line[6]){
      for (6..$#line){
	$extension .= $line[$_]."\t";
      }
      $extension = substr($extension,0,-1);
    }
    my $feat = Bio::ViennaNGS::Feature->new(chromosome=>$chromosome,
					    start=>$start,
					    end=>$end,
					    name=>$name,
					    score=>$score,
					    strand=>$strand,
					    extension=>$extension);
    push @featurelist, $feat;
  }
  return (\@featurelist);
}

sub fetch_chrom_sizes{
  my $species = shift;
  my %sizes;
  my @chromsize;
  my $this_function = (caller(0))[3];

  my $test_fetchChromSizes = can_run('fetchChromSizes') or
    say "ERROR [$this_function] fetchChromSizes utility not found";

  my $cmd = "fetchChromSizes $species";
  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(command => $cmd, verbose => 0);
  if ($success){
    @chromsize = @{$stdout_buf};
  }
  else{
    print STDERR "Using UCSCs fetchChromSizes failed, trying alternative mysql fetch!\n";
    my $test_fetchChromSizes = can_run('mysql') or
      die "ERROR [$this_function] mysql utility not found";
    $cmd = "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from $species.chromInfo\"";  ### Alternative to UCSC fetchChromSizes, has mysql dependency
    my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf )  = run(command => $cmd, verbose => 0);
    if ($success){
      @chromsize = @{$stdout_buf};
    }
    else{
      print STDERR "ERROR [$this_function] External command call unsuccessful\n";
      print STDERR "ERROR: this is what the command printed:\n";
      print join "", @$full_buf;
      croak $!;
      die "Fetching of chromosome sizes failed, please either download fetchChromSizes from the UCSC script collection, or install mysql!\n";
    }
  }

  foreach (@chromsize){
    chomp($_);
    foreach (split(/\n/,$_)){
      my ($chr,$size)=split (/\t/,$_);
      $sizes{$chr}=$size;
    }
  }

  return(\%sizes);
}

1;
__END__


=head1 NAME

Bio::ViennaNGS::Util - Utility routines for Next-Generation Sequencing data
analysis

=head1 SYNOPSIS

  use Bio::ViennaNGS::Util;

  # make bigWig from BED or BAM
  $type = "bam";
  $strand = "+";
  $bwfile = bed_or_bam2bw($type,$infile,$cs_in,$strand,$destdir,$wantnorm,$size_p,$scale,$logfile);

  # make bigBed from BED
  my $bb = bed2bigBed($bed_in,$cs_in,$destdir,$logfile);

  # sort a BED file 
  sortbed($bed_in,$destdir,$bed_out,$rm_orig,$logfile)

  # compute transcript abundance in TPM
  $meanTPM = computeTPM($sample,$readlength);

  # parse a bedtools multicov compatible file
  $conds = parse_multicov($infile);

  # write bedtools multicov compatible file
  write_multicov("TPM",$destdir,$basename);

=head1 DESCRIPTION

Bio::ViennaNGS::Util is a collection of utility subroutines for
building efficient Next-Generation Sequencing (NGS) data analysis
pipelines.

=head2 ROUTINES

=over

=item bed_or_bam2bw($type,$infile,$chromsizes,$strand,$dest,$want_norm,$size,$scale,$log)

Creates stranded, normalized BigWig coverage profiles from
strand-specific BAM or BED files (provided via C<$infile>). The
routine expects a file type 'bam' or 'bed' via the C<$type>
argument. C<$chromsizes> is the chromosome.sizes files, C<$strand> is
either "+" or "-" and C<$dest> contains the output path for
results. For normlization of bigWig profiles, additional attributes
are required: C<$want_norm> triggers normalization with values 0 or
1. C<$size> is the number of fragments/elements in the BAM or BED file
and C<$scale> gives the number to which data is normalized (ie. every
bedGraph entry is multiplied by a factor (C<$scale>/C<$size>). C<$log>
is expected to contain either the full path and file name of log file
or 'undef'. The routine returns the full file name of the newly
generated bigWig file.

While this routine can handle non-straned BAM/BED files (in which case
C<$strand> should be set to "+" and hence all coverage profiles will
be created with a positive sign, even if they map to the negative
strand), usage of strand-specific data is highly recommended. For BAM
file, this is easily achieved by calling the L<bam_split> routine (see
above) prior to this one, thus creating dedicated BAM files containing
exclusively reads mapped to the positive or negative strand,
respectively.

It is important to know that this routine B<does not extract> reads
mapped to either strand from a non-stranded BAM/BED file if the
C<$strand> argument is given. It rather adjusts the sign of B<all>
mapped reads/features in a BAM/BED file and then creates bigWig
files. See the L<split_bam> routine for extracting reads mapped to
either strand.

Stranded bigWigs can easily be visualized via
L<TrackHubs|http://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html>
in the UCSC Genome Browser. Internally, the conversion from BAM/BED to
bigWig is accomplished via two third-party applications:
L<genomeCoverageBed|http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html>
and
L<bedGraphToBigWig|http://hgdownload.cse.ucsc.edu/admin/exe/>. Intermediate
bedGraph files are removed automatically once the bigWig files are
ready.

=item sortbed($infile,$dest,$outfile,$rm_orig,$log)

Sorts BED file C<$infile> with F<bedtools sortt>. C<$dest> and
C<outfile> name path and filename of the resulting sorted BED
file. C<$rm_infile> is either 1 or 0 and indicated whether the
original C<$infile> should be deleted. C<$log> holds path and name of
log file.

=item bed2bigBed($infile,$chromsizes,$dest,$log)

Creates an indexed bigBed file from a BED file. C<$infile> is the BED
file to be transformed, C<$chromsizes> is the chromosome.sizes file
and C<$dest> contains the output path for results. C<$log> is the name
of a log file, or undef if no logging is reuqired. A '.bed', '.bed6'
or '.bed12' suffix in C<$infile> will be replace by '.bb' in the
output. Else, the name of the output bigBed file will be the value of
C<$infile> plus '.bb' appended.

The conversion from BED to bigBed is done by a third-party utility
(bedToBigBed), which is executed by IPC::Cmd.

=item computeTPM($featCount_sample,$rl)

Computes expression in Transcript per Million (TPM) [Wagner
et.al. Theory Biosci. (2012)]. C<$featCount_sample> is a reference to
a Hash of Hashes data straucture where keys are feature names and
values hold a hash that must at least contain length and raw read
counts. Practically, C<$featCount_sample> is represented by _one_
element of C<@featCount>, which is populated from a multicov file by
C<parse_multicov()>. C<$rl> is the read length of the sequencing run.

Returns the mean TPM of the processed sample, which is invariant among
samples. (TPM models relative molar concentration and thus fulfills
the invariant average criterion.)

=item parse_multicov($file)

Parse a bedtools multicov (multiBamCov) file, i.e. an extended BED6
file, into an Array of Hash of Hashes data structure
(C<@featCount>). C<$file> is the input file. Returns the number of
samples present in the multicov file, ie. the numner of columns
extending the canonical BED6 columns in the input multicov file.

=item write_multicov($item,$dest,$base_name)

Write C<@featCount> data to a bedtools multicov (multiBamCov)-type
file. C<$item> specifies the type of information from C<@featCount>
HoH entries, e.g. TPM or RPKM. These values must have been computed
and inserted into C<@featCount> beforehand by
e.g. C<computeTPM()>. C<$dest> gives the absolute path and
C<$base_name> the basename (will be extended by C<$item>.csv) of the
output file.

=back

=head1 DEPENDENCIES

=over 7

=item  L<Bio::Perl> >= 1.00690001

=item  L<BIO::DB::Sam> >= 1.39

=item  L<File::Basename>

=item  L<File::Temp>

=item  L<Path::Class>

=item  L<IPC::Cmd>

=item  L<Carp>

=back

L<Bio::ViennaNGS> uses third-party tools for computing intersections
of BED files: F<bedtools intersect> from the
L<BEDtools|http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html>
suite is used to compute overlaps and F<bedtools sort> is used to sort
BED output files. Make sure that those third-party utilities are
available on your system, and that hey can be found and executed by
the Perl interpreter. We recommend installing the latest version of
L<BEDtools|https://github.com/arq5x/bedtools2> on your system.

=head1 AUTHORS

=over

=item Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=item Jörg Fallmann E<lt>fall@tbi.univie.ac.atE<gt>

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
