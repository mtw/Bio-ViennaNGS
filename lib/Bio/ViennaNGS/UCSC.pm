# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-26 17:28:26 mtw>

package Bio::ViennaNGS::UCSC;

use Bio::ViennaNGS;
use Exporter;
use strict;
use warnings;
use Template;
use Cwd;
use File::Basename;
use File::Find;
use IPC::Cmd qw(can_run run);
use File::Share ':all';
use Path::Class;
use Data::Dumper;
use Carp;
use Bio::ViennaNGS::Util qw(fetch_chrom_sizes bed2bigBed);
use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");

our @ISA = qw(Exporter);

our @EXPORT_OK = qw( make_assembly_hub make_track_hub retrieve_color
make_track make_multi_bigwig_container_track
make_bigwig_container_track retrieve_bigwig_tracks
retrieve_bigbed_url_tracks parse_fasta_header valid_ncbi_accession
modify_fasta_header make_group retrieve_chromosome_size
write_chromosome_size_file convert_tracks);

our @EXPORT = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub make_assembly_hub{
  my ($fasta_path, $filesdir, $basedir, $baseURL, $big_wig_ids, $log) = @_;
  my ($basename,$dir,$ext);
  my $this_function = (caller(0))[3];
  #check arguments
  croak ("ERROR [$this_function] \$fasta_path does not exist\n") 
    unless (-e $fasta_path);
  croak ("ERROR [$this_function] \$basedir does not exist\n") 
    unless (-d $basedir);
  croak ("ERROR [$this_function]: no URL (network location for upload to UCSC) provided") 
    unless(defined $baseURL);

  unless ($baseURL =~ /\/$/) { $baseURL .= "/"; }

  my $tmp_path = dist_file('Bio-ViennaNGS', "hub.txt" );
  ($basename,$dir,$ext) = fileparse($tmp_path,qr/\..*/);
  my $template_path = dir($dir,"template-AssemblyHub");

  croak ("ERROR [$this_function] template directory not found\n") 
    unless (-d $template_path);
  my $faToTwoBit = can_run('faToTwoBit') or
    croak ("ERROR [$this_function] faToTwoBit is not installed!");

  # bedfiles path
  my $parsedHeaderRef = parse_fasta_header($fasta_path,$this_function);
  my @parsedHeader = @$parsedHeaderRef;
  my $unchecked_accession = $parsedHeader[0];
  my $scientificName = $parsedHeader[1];
  my $accession = valid_ncbi_accession($unchecked_accession);
  # create assembly hub directory structure
  my $assembly_hub_name = "assemblyHub";
  my $assembly_hub_directory = dir($basedir, $assembly_hub_name);
  my $genome_assembly_name = $accession;
  my $genome_assembly_directory = dir($assembly_hub_directory,$genome_assembly_name);
  mkdir $assembly_hub_directory;
  mkdir $genome_assembly_directory;

  if (defined $log){
    open(LOG, ">>", $log) or croak "$!";
    print LOG "LOG Base directory:          $assembly_hub_directory\n";
    print LOG "LOG Assembly Hub directory:  $genome_assembly_directory\n";
    close(LOG);
  }

  #2-bit fasta file conversion
  my $fa_modified = file($assembly_hub_directory, $accession.".fa");
  modify_fasta_header($fasta_path,$fa_modified,$accession);
  my $twoBit = file($genome_assembly_directory, $accession.".2bit");
  my $fastaToTwobit_cmd = "$faToTwoBit $fa_modified $twoBit";

  if (defined $log){
    open(LOG, ">>", $log) or croak "$!";
    print LOG "LOG [$this_function] $fastaToTwobit_cmd\n";
    close(LOG);
  }

  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $fastaToTwobit_cmd, verbose => 0 );
  if( !$success ) {
    print STDERR "ERROR [$this_function] External command call unsuccessful\n";
    print STDERR "ERROR: this is what the command printed:\n";
    print join "", @$full_buf;
    croak $!;
  }

  #template definition
  my $template = Template->new({
                                INCLUDE_PATH => ["$template_path"],
                                RELATIVE=>1,
  });

  #construct hub.txt
  my $hubtxt_path = file($assembly_hub_directory,"hub.txt")->stringify;
  my $hubtxt_file = "hub.txt";
  my $hubtxt_vars =
    {
     hubName => $accession,
     shortLabel => $accession,
     longLabel => $accession,
     genomesFile => "genome.txt",
     email => 'email',
     descriptionURL => "$baseURL" . "description.html"
    };
  $template->process($hubtxt_file,$hubtxt_vars,"$hubtxt_path") ||
    croak "Template process failed: ", $template->error(), "\n";

  #construct genome.txt
  my $genometxt_path = file($assembly_hub_directory, "genome.txt")->stringify;
  my $genometxt_file = "genome.txt";
  my $genometxt_vars =
    {
     genome => $accession,
     trackDb => file($accession, "trackDb.txt"),
     groups => file($accession, "groups.txt"),
     description => "$accession",
     twoBitPath => file($accession,$accession.".2bit"),
     organism => "organism",
     defaultPos => $accession,
     orderKey => "10",
     scientificName => "$scientificName",
     htmlPath => file($accession,"description.html")
    };
  $template->process($genometxt_file,$genometxt_vars,$genometxt_path) or
    croak "Template process failed: ", $template->error(), "\n";

  #construct description.html
  my $description_html_path = file($genome_assembly_directory, "description.html")->stringify;
  my $description_html_file = "description.html";
  my $description_html_vars =
    {
     imageLink  => "imageLink",
     imageSource => "imageSource",
     imageAlternative => "imageAlternative",
     taxonomicName => "taxonomicName",
     imageOrigin => "imageOrigin",
     imageOriginDescription => "imageOriginDescription",
     ucscId => "ucscId",
     sequencingId => "sequencingId",
     assemblyDate => "assemblyDate",
     genbankAccessionId => "genbankAccessionId",
     ncbiGenomeInformationLink => "ncbiGenomeInformationLink",
     ncbiGenomeInformationDescription => "ncbiGenomeInformationDescription",
     ncbiAssemblyInformationLink => "ncbiAssemblyInformationLink",
     ncbiAssemblyInformationDescription => "ncbiAssemblyInformationDescription",
     bioProjectInformationLink => "bioProjectInformationLink",
     bioProjectInformationDescription => "bioProjectInformationDescription",
     sequenceAnnotationLink => "sequenceAnnotationLink"
    };
  $template->process($description_html_file,$description_html_vars,$description_html_path) or
    croak "Template process failed: ", $template->error(), "\n";

  my $groups = make_group("annotation", "Annotation", "1", "0");

  #construct group.txt
  my $group_txt_path = file($genome_assembly_directory, "groups.txt")->stringify;
  my $group_txt_file = 'groups.txt';
  my $group_txt_vars = 
    {
     groups  => "$groups",
    };
  $template->process($group_txt_file,$group_txt_vars,$group_txt_path) or
    croak "Template process failed: ", $template->error(), "\n";

  #construct big bed
  my $tracksList;
  #Bigbeds are only created from infolder
  unless($filesdir =~ /-/){
    my $chromosome_size = retrieve_chromosome_size($fasta_path);
    my $chromosome_size_filepath = file($genome_assembly_directory,"$accession.chrom.sizes");
    write_chromosome_size_file($chromosome_size_filepath,$accession,$chromosome_size);
    convert_tracks($filesdir, $genome_assembly_directory, $chromosome_size_filepath, $log);
    my @trackfiles = retrieve_tracks($genome_assembly_directory, $baseURL, $assembly_hub_name, $accession);

    foreach my $track (@trackfiles){
      my $trackString = make_track(@$track);
      $tracksList .= $trackString;
    }
  }

  #big wigs
  my $bigwig_tracks_string = "";
  unless($big_wig_ids=~/^-$/){
    $bigwig_tracks_string = retrieve_bigwig_tracks($genome_assembly_directory, $baseURL, $assembly_hub_name, $accession, $big_wig_ids);
  }
  $tracksList .= $bigwig_tracks_string;
  #construct trackDb.txt
  my $trackDb_txt_path = file($genome_assembly_directory, "trackDb.txt")->stringify;
  my $trackDb_txt_file = 'trackDb.txt';
  my $trackDb_txt_vars =
    {
     tracks  => "$tracksList"
    };
  $template->process($trackDb_txt_file,$trackDb_txt_vars,$trackDb_txt_path) or
    croak "Template process failed: ", $template->error(), "\n";

  if (defined $log){
    open(LOG, ">>", $log) or croak "$!";
    print LOG "LOG Assembly Hub created\n";
    close(LOG);
  }
}

sub make_track_hub{
  my ($species, $basedir, $baseURL, $big_bed_ids, $big_wig_ids, $log) = @_;
  my ($basename,$dir,$ext);
  my $this_function = (caller(0))[3];

  #check arguments
  croak ("ERROR [$this_function] \no species provided\n")
    unless ($species);
  croak ("ERROR [$this_function] \$basedir does not exist\n")
    unless (-d $basedir);
  croak ("ERROR [$this_function]: no URL (network location for upload to UCSC) provided")
    unless(defined $baseURL);

  if (defined $log){
    open(LOG, ">>", $log) or croak "$!";
  }

  unless ($baseURL =~ /\/$/) { $baseURL .= "/"; }

  my $tmp_path = dist_file('Bio-ViennaNGS', "hub.txt" );
  ($basename,$dir,$ext) = fileparse($tmp_path,qr/\..*/);
  my $template_path = dir($dir,"template-TrackHub");

  croak ("ERROR [$this_function] template directory not found\n") 
    unless (-d $template_path);
  my $faToTwoBit = can_run('faToTwoBit') or
    croak ("ERROR [$this_function] faToTwoBit is not installed!");

  # create track hub directory structure
  my $track_hub_name = "trackHub";
  my $track_hub_directory = dir($basedir, $track_hub_name);
  my $genome_assembly_name = $species;
  my $genome_assembly_directory = dir($track_hub_directory,$genome_assembly_name);
  mkdir $track_hub_directory;
  mkdir $genome_assembly_directory;
  if (defined $log){
    print LOG "LOG Base directory:          $track_hub_directory\n";
    print LOG "LOG Track Hub directory:  $genome_assembly_directory\n";
  }

  #template definition
  my $template = Template->new({
                                INCLUDE_PATH => ["$template_path"],
                                RELATIVE=>1,
  });

  #construct hub.txt
  my $hubtxt_path = file($track_hub_directory,"hub.txt")->stringify;
  my $hubtxt_file = "hub.txt";
  my $hubtxt_vars =
    {
     hubName => $species,
     shortLabel => $species,
     longLabel => $species,
     genomesFile => "genome.txt",
     email => 'email',
     descriptionURL => "$baseURL" . "description.html"
    };
  $template->process($hubtxt_file,$hubtxt_vars,"$hubtxt_path") ||
    croak "Template process failed: ", $template->error(), "\n";

  #construct genome.txt
  my $genometxt_path = file($track_hub_directory, "genome.txt")->stringify;
  my $genometxt_file = "genome.txt";
  my $genometxt_vars =
    {
     genome => $species,
     trackDb => file($species, "trackDb.txt"),
     groups => file($species, "groups.txt"),
     description => "$species",
     twoBitPath => file($species,$species.".2bit"),
     organism => "organism",
     defaultPos => $species,
     orderKey => "10",
     scientificName => "scientificName",
     htmlPath => file($species,"description.html")
    };
  $template->process($genometxt_file,$genometxt_vars,$genometxt_path) or
    croak "Template process failed: ", $template->error(), "\n";

  my $tracksList;
  #bigbeds are only created from urls
  my $bigbeds_tracks_string = "";
  unless($big_bed_ids=~/^-$/){
    $bigbeds_tracks_string = retrieve_bigbed_url_tracks($genome_assembly_directory, $baseURL, $track_hub_name, $big_bed_ids);
  }
  $tracksList .= $bigbeds_tracks_string;
  #big wigs
  my $bigwig_tracks_string = "";
  unless($big_wig_ids=~/^-$/){
    $bigwig_tracks_string = retrieve_bigwig_tracks($genome_assembly_directory, $baseURL, $track_hub_name, $species, $big_wig_ids);
  }
  $tracksList .= $bigwig_tracks_string;
  #construct trackDb.txt
  my $trackDb_txt_path = file($genome_assembly_directory, "trackDb.txt")->stringify;
  my $trackDb_txt_file = 'trackDb.txt';
  my $trackDb_txt_vars =
    {
     tracks  => "$tracksList"
    };
  $template->process($trackDb_txt_file,$trackDb_txt_vars,$trackDb_txt_path) or
    croak "Template process failed: ", $template->error(), "\n";

  if (defined $log){
    print LOG "LOG Track Hub created\n";
    close(LOG);
  }
}

sub convert_tracks{
  my ($filesdir,$genome_assembly_directory,$chromosome_size_filepath,$log) = @_;
  my $this_function = (caller(0))[3];
  my @bedfiles = ();

  find( sub {
	  return unless -f;
	  return unless /\.bed$/;
	  push @bedfiles, $File::Find::name;
	} , $filesdir);

  if (defined $log){
    open(LOG, ">>", $log) or croak $!;
    foreach my $bedfile (@bedfiles){
      print LOG "LOG [$this_function] found file $bedfile\n";
    }
    close(LOG);
  }

  foreach my $bedfile (@bedfiles){
    bed2bigBed($bedfile,$chromosome_size_filepath,$genome_assembly_directory,$log);
  }
  return 1;
}

sub retrieve_tracks{
  my ($directoryPath,$baseURL,$assembly_hub_name,$accession) = @_;
  my $currentDirectory = getcwd;
  chdir $directoryPath or croak $!;
  my @trackfiles = <*.bb>;
  my @tracks;
  my $counter = 0;
  foreach my $trackfile (@trackfiles){
    my $filename = $trackfile;
    $filename =~ s/.bb$//;
    my @filenameSplit = split(/\./, $filename);
    my $id = lc($filename);
    my $tag = $id;
    my $track = $id . "_bed";
    my $dir = dir($baseURL,$assembly_hub_name,$accession);
    my $bigDataUrl = file($trackfile);
    my $shortLabel = $id;
    my $longLabel = $id;
    my $type = "bigBed 12 .";
    my $autoScale = "off";
    my $bedNameLabel = "Gene Id";
    my $searchIndex = "name";
    my $colorByStrand = retrieve_color($counter);
    my $visibility = "pack";
    my $group = "annotation";
    my $priority = "10";
    my @track = ($tag,$track,$bigDataUrl,$shortLabel,$longLabel,$type,$autoScale,$bedNameLabel,$searchIndex,$colorByStrand,$visibility,$group,$priority);
    my $trackreference = \@track;
    push(@tracks, $trackreference);
    $counter++;
  }
  chdir $currentDirectory or croak $!;
  return @tracks;
}

sub retrieve_bigbed_url_tracks{
  my ($directoryPath,$baseURL,$assembly_hub_name,$bigbedpaths) = @_;
  my $currentDirectory = getcwd;
  chdir $directoryPath or croak $!;
  my @big_bed_array = split('#', $bigbedpaths);
  my $bigbedtracks = "";
  my $counter = 0;
  foreach my $bigbed_entry (@big_bed_array){
      #construct single big bed track
      my ($basename,$dir,$ext) = fileparse($bigbed_entry,qr/\..*/);
      my $id = lc($basename);
      my $tag = $id . $ext . "_bed";
      my $track = $id . $ext  ."_bed";
      my $bigDataUrl = $bigbed_entry;
      my $shortLabel = $id;
      my $longLabel = $id;
      my $type = "bigBed 12 .";
      my $autoScale = "off";
      my $bedNameLabel = "Gene Id";
      my $searchIndex = "name";
      my $colorByStrand = retrieve_color($counter);
      my $visibility = "pack";
      my $group = "annotation";
      my $priority = "10";
      my $track_string = make_track($tag, $track, $bigDataUrl, $shortLabel, $longLabel, $type, $autoScale, $bedNameLabel, $searchIndex, $colorByStrand, $visibility, $group, $priority);
      $bigbedtracks .= $track_string;
      $counter++;
  }
  chdir $currentDirectory or croak $!;
  return $bigbedtracks;
}

sub retrieve_bigwig_tracks{
  my ($directoryPath,$baseURL,$assembly_hub_name,$accession,$bigwigpaths) = @_;
  my $currentDirectory = getcwd;
  chdir $directoryPath or croak $!;
  my @big_wig_array = split('\#', $bigwigpaths);
  my $bigwigtracks = "";
  foreach my $bigwig_entry (@big_wig_array){
    if($bigwig_entry=~/,/){
      #multi bigwig container detected
      my @big_wig_pair = split(',', $bigwig_entry);
      my $pos_wig=$big_wig_pair[0];
      my $neg_wig=$big_wig_pair[1];
      my ($basename1,$dir1,$ext1) = fileparse($pos_wig,qr/\..*/);
      my ($basename2,$dir2,$ext2) = fileparse($pos_wig,qr/\..*/);
      #construct container
      my $id = lc($basename1);
      my $tag = $id . "_bw";
      my $track = $id . "_bw";
      my $shortLabel = $tag;
      my $longLabel = $tag;
      my $type = "bigWig";
      my $autoScale = "on";
      my $visibility = "full";
      my $priority = "1500";
      my $container_string = make_multi_bigwig_container_track($tag, $track, $shortLabel, $longLabel, $type, $autoScale, $visibility, $priority);
      $bigwigtracks .= $container_string;
      #construct positive track
      my $track1 = $track . "_pos";
      my $bigDataUrl1 = $pos_wig;
      my $shortLabel1 = $track1;
      my $longLabel1 = $track1;
      my $type1 = "bigWig";
      my $parent1 = $track;
      my $color1 = "133,154,0";
      my $track1_string = make_bigwig_container_track($track1, $bigDataUrl1, $shortLabel1, $longLabel1, $type1, $parent1, $color1);
      $bigwigtracks .= $track1_string;
      #construct negative track
      my $track2 = $track . "_neg";
      my $bigDataUrl2 = $neg_wig;
      my $shortLabel2 = $track2;
      my $longLabel2 = $track2;
      my $type2 = "bigWig";
      my $parent2 = $track;
      my $color2 = "220,51,47";
      my $track2_string = make_bigwig_container_track($track2, $bigDataUrl2, $shortLabel2, $longLabel2, $type2, $parent2, $color2);
      $bigwigtracks .= $track2_string;
    }else{
      #construct single big wig track
      my ($basename,$dir,$ext) = fileparse($bigwig_entry,qr/\..*/);
      my $id = lc($basename);
      my $tag = $id . "_bw";
      my $track = $id . "_bw";
      my $bigDataUrl = $bigwig_entry;
      my $shortLabel = $tag;
      my $longLabel = $tag;
      my $type = "bigWig";
      my $autoScale = "on";
      my $visibility = "full";
      my $priority = "1500";
      my $color = "38,140,210";
      my $track_string = make_bigwig_track($tag, $track, $bigDataUrl, $shortLabel, $longLabel, $type, $autoScale, $visibility, $priority, $color);
      $bigwigtracks .= $track_string;
    }
  }
  chdir $currentDirectory or croak $!;
  return $bigwigtracks;
}


sub retrieve_color{
  my $counter = shift;
  my $digitnumber = length($counter);
  my $colorcode;
  if($digitnumber>1){
    $colorcode = $counter % 10;
  }else{
    $colorcode = $counter;
  }
  my @color;
  #green
  $color[0] = "133,154,0 133,154,0";
  #cyan
  $color[1] = "42,162,152 42,162,152";
  #blue
  $color[2] = "38,140,210 38,140,210";
  #violet
  $color[3] = "108,114,196 108,114,196";
  #magenta
  $color[4] = "211,55,130 211,55,130";
  #red
  $color[5] = "220,51,47 220,51,47";
  #orange
  $color[6] = "203,76,22 203,76,22";
  #yellow
  $color[7] = "181,138,0 181,138,0";
  #black
  $color[8] = "0,44,54 0,44,54";
  #grey
  $color[9] = "88,111,117 88,111,117";

  return $color[$colorcode];
}

sub make_group{
  my ($name, $label, $priority, $defaultIsClosed) = @_;
  my $group ="name $name\nlabel $label\npriority $priority\ndefaultIsClosed $defaultIsClosed\n";
  return $group;
}

sub make_track{
  my ($tag, $track, $bigDataUrl, $shortLabel, $longLabel, $type, $autoScale, $bedNameLabel, $searchIndex, $colorByStrand, $visibility, $group, $priority) = @_;
  my $trackEntry ="#$tag\ntrack $track\nbigDataUrl $bigDataUrl\nshortLabel $shortLabel\nlongLabel $longLabel\ntype $type\nautoScale $autoScale\nbedNameLabel $bedNameLabel\nsearchIndex $searchIndex\ncolorByStrand $colorByStrand\nvisibility $visibility\ngroup $group\npriority $priority\n\n";
  return $trackEntry;
}

sub make_multi_bigwig_container_track{
  my ($tag, $track, $shortLabel, $longLabel, $type, $autoScale, $visibility, $priority) = @_;
  my $trackEntry ="#$tag\ntrack $track\ncontainer multiWig\nnoInherit on\nshortLabel $shortLabel\nlongLabel $longLabel\ntype $type\nconfigureable on\nvisibility $visibility\naggregate transparentOverlay\nshowSubtrackColorOnUi on\nautoScale $autoScale\nwindowingFunction maximum\npriority $priority\nalwaysZero on\nyLineMark 0\nyLineOnOff on\nmaxHeightPixels 125:125:11\n\n";
  return $trackEntry;
}

sub make_bigwig_container_track{
  my ($track, $bigDataUrl, $shortLabel, $longLabel, $type, $parent, $color) = @_;
  my $trackEntry = "track $track\nbigDataUrl $bigDataUrl\nshortLabel $shortLabel\nlongLabel $longLabel\ntype $type\nparent $parent\ncolor $color\n\n";
  return $trackEntry;
}

sub make_bigwig_track{
  my ($tag, $track, $bigDataUrl, $shortLabel, $longLabel, $type, $autoScale, $visibility, $priority, $color) = @_;
  my $trackEntry ="#$tag\ntrack $track\nbigDataUrl $bigDataUrl\nshortLabel $shortLabel\nlongLabel $longLabel\ntype $type\nvisibility $visibility\nautoScale $autoScale\npriority $priority\nalwaysZero on\nyLineMark 0\nyLineOnOff on\nmaxHeightPixels 125:125:11\ncolor $color\n\n";
  return $trackEntry;
}

sub valid_ncbi_accession{
  # receives a NCBI accession ID, with or without version number
  # checks for validity and returns the accession number as is
  my $acc = shift;
  if ($acc =~ /^(N[CSTZ]\_[A-Z]*\d{6})\.\d+?$/){
    return $acc;
  }
  elsif ($acc =~ /^(N[CSTZ]\_[A-Z]*\d{6})$/){
    return $1;
  }
  else {
    return 0;
  }
}

sub parse_fasta_header{
  my $filepath = shift;
  my $this_function = (caller(0))[3];
  open my $file, '<', "$filepath" or die $!;
  my $fastaheader = <$file>;
  chomp $fastaheader;
  close $file;
  my @ids = ();
  #>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655
  if($fastaheader=~/^>gi/){
    my @headerfields = split(/\|/, $fastaheader);
    my $accession = $headerfields[3];
    my $scientificName = $headerfields[4];
    push(@ids,$accession);
    push(@ids,$scientificName);
    return \@ids;
  }else{
    $fastaheader=~s/^>//;
    carp "INFO [$this_function] It looks like your input Fasta header does not contain a valid NCBI accession number. Continuing with:\n $fastaheader\n"
      unless ( valid_ncbi_accession($fastaheader) );
    push(@ids,$fastaheader);
    push(@ids,"scientific name not set");
    return \@ids;
  }
}

sub write_chromosome_size_file{
  my $filepath = shift;
  my $chromosome_name = shift;
  my $chromosome_size = shift;
  my $entry = $chromosome_name . "\t" . $chromosome_size . "\n";
  open CHROMFILE, '>', "$filepath" or die $!;
  print CHROMFILE $entry;
  close CHROMFILE;
  return 1;
}

sub write_chromosome_sizes_file{
  my $filepath = shift;
  my $chromosome_sizes_reference = shift;
  my %chromosome_sizes = %{$chromosome_sizes_reference};
  open CHROMFILE, '>', "$filepath" or die $!;
  foreach my $chromosome_name ( keys %chromosome_sizes){
    my $chromosome_size = $chromosome_sizes{$chromosome_name};
    my $entry = $chromosome_name . "\t" . $chromosome_size . "\n";
    print CHROMFILE $entry;
  }
  close CHROMFILE;
  return 1;
}

sub retrieve_chromosome_size{
  my $inputFilepath = shift;
  open INFILE, '<', $inputFilepath;
  my @newfasta;
  my $chromosome_size = 0;
  my $header_skipped = 0;
  while (<INFILE>) {
    if($header_skipped){
      chomp;
      $chromosome_size += length($_);
    }else{
      $header_skipped = 1;
    }
  }
  close INFILE;
  return $chromosome_size;
}

sub modify_fasta_header{
  my $inputFilepath = shift;
  my $outputFilepath = shift;
  my $header = shift;

  open INFILE, '<', "$inputFilepath";
  my @newfasta;
  while (<INFILE>) {
   push(@newfasta, $_);
  }
  close INFILE;
  #@newfasta  = @newfasta[ 1 .. $#newfasta ];
  shift @newfasta;
  unshift(@newfasta,">".$header."\n");

  open OUTFILE, '>', "$outputFilepath";
  foreach my $line (@newfasta) {
   print OUTFILE $line;
  }
  close OUTFILE;
  return 1;
}

1;
__END__

=head1 NAME

Bio::ViennaNGS::UCSC - Perl extension for easy UCSC Genome Browser
integration.

=head1 SYNOPSIS

  use Bio::ViennaNGS::UCSC;

=head1 DESCRIPTION

L<Bio::ViennaNGS::UCSC> facilitates routine tasks for automatic
visualization of next-generation sequencing data in the UCSC Genome
Browser.

=head2 EXPORT

Routines: 
  make_assembly_hub, make_track_hub

Variables:
  none

=head3 make_assembly_hub()

Build an UCSC Genome Browser Assembly Hub. This function takes 4
parameters:

=over

=item 1 absolute path of input fasta file(e.g. /home/user/input.fa)

=item 2 path to the ouput directory (e.g. /home/user/assemblyhubs/)

=item 3 base URL under which the new Assembly Hub will be available on
the Internet, e.g. http://www.foo.com/folder/. The generated Assembly
Hub files must be copied to the corresponding file system location on
the Web Server and the URL readable for the UCSC Genome Browser server.

=item 4 path for the log file (/home/user/logs/assemblyhubconstructionlog)

=back

=head3 make_track_hub()

Build track hubs for the UCSC genome browser.
This function takes 4 parameters:

=over

=item 1 chromosome id as used in existing ucsc assembly hub (e.g. chr1)

=item 2 path to the ouput directory (e.g. /home/user/assemblyhubs/)

=item 3 base URL where the output folder will be placed for upload to the UCSC genome browser (e.g. http://www.foo.com/folder/)

=item 4 path for the log file (/home/user/logs/assemblyhubconstructionlog)

=back

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=back

=head1 AUTHORS

=over

=item Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=item Florian Eggenhofer, E<lt>florian.eggenhofer@univie.ac.atE<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2017 Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


=cut
