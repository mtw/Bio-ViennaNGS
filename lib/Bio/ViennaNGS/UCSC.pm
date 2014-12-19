# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:34:07 mtw>

package Bio::ViennaNGS::UCSC;

use Exporter;
use version; our $VERSION = qv('0.12_07');
use strict;
use warnings;
use Template;
use Cwd;
use File::Basename;
use IPC::Cmd qw(can_run run);
use File::Share ':all';
use Path::Class;
use Data::Dumper;
use Carp;
use Bio::ViennaNGS::UCSC;

our @ISA = qw(Exporter);

our @EXPORT_OK = qw( make_assembly_hub make_track_hub );

our @EXPORT = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub make_assembly_hub{
  my ($fasta_path, $filesdir, $basedir, $baseURL, $log) = @_;
  my ($basename,$dir,$ext);
  my $this_function = (caller(0))[3];

  #check arguments
  croak ("ERROR [$this_function] \$fasta_path does not exist\n") 
    unless (-e $fasta_path);
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
  my $template_path = dir($dir,"template-AssemblyHub");

  croak ("ERROR [$this_function] template directory not found\n") 
    unless (-d $template_path);
  my $faToTwoBit = can_run('faToTwoBit') or
    croak ("ERROR [$this_function] faToTwoBit is not installed!");

  my $bedToBigBed = can_run('bedToBigBed') or
    croak ("ERROR [$this_function] bedToBigBed is not installed!");

  # bedfiles path
  my @parsedHeader = parse_fasta_header($fasta_path);
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
    print LOG "LOG Base directory:          $assembly_hub_directory\n";
    print LOG "LOG Assembly Hub directory:  $genome_assembly_directory\n";
  }

  #2-bit fasta file conversion
  my $fa_modified = file($assembly_hub_directory, $accession.".fa");
  modify_fasta_header($fasta_path,$fa_modified,$accession);
  my $twoBit = file($genome_assembly_directory, $accession.".2bit");
  my $fastaToTwobit_cmd = "$faToTwoBit $fa_modified $twoBit";
  if (defined $log){  print LOG "LOG [$this_function] $fastaToTwobit_cmd\n";}

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


  my $chromosome_size = retrieve_chromosome_size($fasta_path);
  my $chromosome_size_filepath = file($genome_assembly_directory,"$accession.chrom.sizes");
  write_chromosome_size_file($chromosome_size_filepath,$accession,$chromosome_size);
  convert_tracks($filesdir, $genome_assembly_directory, $accession, $bedToBigBed, $chromosome_size_filepath);
  my @trackfiles = retrieve_tracks($genome_assembly_directory, $baseURL, $assembly_hub_name, $accession);

  my $tracksList;
  foreach my $track (@trackfiles){
    my $trackString = make_track(@$track);
    $tracksList .= $trackString;
  }

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
    print LOG "LOG Assembly Hub created\n";
    close(LOG);
  }
}

sub make_track_hub{
  my ($species, $filesdir, $basedir, $baseURL, $chrom_sizes_file, $chrom_size_file, $log) = @_;
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

  my $bedToBigBed = can_run('bedToBigBed') or
    croak ("ERROR [$this_function] bedToBigBed is not installed!");

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

  if(-e $chrom_sizes_file){
    convert_tracks($filesdir, $genome_assembly_directory, $species, $bedToBigBed, $chrom_sizes_file);
  }else{
    my $chromosome_sizes = fetch_chrom_sizes($species);
    my $chromosome_size_filepath = file($genome_assembly_directory,"$species.chrom.sizes");
    write_chromosome_sizes_file($chromosome_size_filepath,$chromosome_sizes);
    convert_tracks($filesdir, $genome_assembly_directory, $species, $bedToBigBed, $chromosome_size_filepath);
  }
  my @trackfiles = retrieve_tracks($genome_assembly_directory, $baseURL, $track_hub_name, $species);

  my $tracksList;
  foreach my $track (@trackfiles){
    my $trackString = make_track(@$track);
    $tracksList .= $trackString;
  }

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
  my ($filesdir,$genome_assembly_directory,$accession,$bedToBigBed,$chromosome_size_filepath) = @_;
  my $currentDirectory = getcwd;
  chdir $filesdir or croak $!;
  my @bedfiles = <*.bed>;
  foreach my $bedfile (@bedfiles){
    my $bed_file_path = file($filesdir,$bedfile);
    my $filename = $bedfile;
    $filename =~ s/.bed$//;
    my $bigbed_file = $filename . ".bb";
    my $bigbed_file_path = file($genome_assembly_directory,$bigbed_file);
    `$bedToBigBed $bed_file_path $chromosome_size_filepath $bigbed_file_path`;
  }
  chdir $currentDirectory or croak $!;
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
    my $fileaccession = $filenameSplit[0];
    my $tag = $filenameSplit[1];
    my $id = lc($tag);
    my $track = "refseq_" . $id;
    my $dir = dir($baseURL,$assembly_hub_name,$accession);
    #my $bigDataUrl = file($dir, $trackfile);
    my $bigDataUrl = file($trackfile);
    my $shortLabel = "RefSeq " . $tag;
    my $longLabel = "RefSeq " . $tag;
    my $type = "bigBed 12 .";
    my $autoScale = "off";
    my $bedNameLabel = "Gene Id";
    my $searchIndex = "name";
    #my $colorByStrand = "100,205,255 55,155,205";
    my $colorByStrand = retrieve_color($counter );
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

sub valid_ncbi_accession{
  # receives a NCBI accession ID, with or without version number
  # returns NCBI accession ID without version number
  my $acc = shift;
  if ($acc =~ /^(N[CST]\_\d{6})\.\d+?$/){
    return $acc; # NCBI accession ID without version
  }
  elsif ($acc =~ /^(N[CST]\_\d{6})$/){
    return $1; # NCBI accession ID without version
  }
  else {
    return 0;
  }
}

sub parse_fasta_header{
  my $filepath = shift;
  open my $file, '<', "$filepath";
  my $fastaheader = <$file>;
  chomp $fastaheader;
  close $file;
  #>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655
  my @headerfields = split(/\|/, $fastaheader);
  my $accession = $headerfields[3];
  my $scientificName = $headerfields[4];
  my @ids;
  push(@ids,$accession);
  push(@ids,$scientificName);
  return @ids;
}

sub write_chromosome_size_file{
  my $filepath = shift;
  my $chromosome_name = shift;
  my $chromosome_size = shift;
  my $entry = $chromosome_name . "\t" . $chromosome_size . "\n";
  open CHROMFILE, '>', "$filepath";
  print CHROMFILE $entry;
  close CHROMFILE;
  return 1;
}

sub write_chromosome_sizes_file{
  my $filepath = shift;
  my $chromosome_sizes_reference = shift;
  my %chromosome_sizes = %{$chromosome_sizes_reference};
  open CHROMFILE, '>', "$filepath";
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

ViennaNGS::UCSC is a Perl extension for managing routine tasks with the
UCSC Genome Browser. It comes with a set of utilities that serve as
reference implementation of the routines implemented in the library. All
utilities are located in the 'scripts' folder.
The main functionality is provided by the make_assembly_hub function.

=head2 EXPORT

Routines: 
  make_assembly_hub, make_track_hub

Variables:
  none

=head3 make_assembly_hub()

Build assembly hubs for the UCSC genome browser.
This function takes 4 parameters:
1. absolute path of input fasta file(e.g. /home/user/input.fa)
2. path to the ouput directory (e.g. /home/user/assemblyhubs/)
3. base URL where the output folder will be placed for upload to the UCSC genome browser (e.g. http://www.foo.com/folder/)
4. path for the log file (/home/user/logs/assemblyhubconstructionlog)

=head3 make_track_hub()

Build track hubs for the UCSC genome browser.
This function takes 4 parameters:
1. chromosome id as used in existing ucsc assembly hub (e.g. chr1)
2. path to the ouput directory (e.g. /home/user/assemblyhubs/)
3. base URL where the output folder will be placed for upload to the UCSC genome browser (e.g. http://www.foo.com/folder/)
4. path for the log file (/home/user/logs/assemblyhubconstructionlog)


=head1 SEE ALSO

perldoc ViennaNGS
perldoc ViennaNGS::AnnoC

=head1 AUTHORS

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>
Florian Eggenhofer, E<lt>florian.eggenhofer@univie.ac.atE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


=cut
