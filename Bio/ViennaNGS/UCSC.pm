# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-22 17:01:47 mtw>

package Bio::ViennaNGS::UCSC;

use Exporter;
use version; our $VERSION = qv('0.01');
use strict;
use warnings;
use Template;
use Cwd;
use File::Basename;
use IPC::Cmd qw[can_run run run_forked];

our @ISA = qw(Exporter);

our @EXPORT_OK = qw( make_assembly_hub  );

our @EXPORT = ();

sub make_assembly_hub{
  my ($fasta_file_path, $assembly_hub_destination_path, $base_URL, $log_path) = @_;
  #check arguments
  die("ERROR [Bio::ViennaNGS::UCSC] \$fasta_file_path does not exist\n") unless (-e $fasta_file_path);
  die("ERROR [Bio::ViennaNGS::UCSC] \$assembly_hub_destination_path does not exist\n") unless (-d $assembly_hub_destination_path);
  die "ERROR [Bio::ViennaNGS::UCSC]: no URL (network location for upload to UCSC) provided" unless(defined $base_URL);
  die("ERROR [Bio::ViennaNGS::UCSC] \$log_path does not exist\n") unless (-e $log_path);
  #ensure that base_URL ends with slash
  $base_URL =~ s!/*$!/!;
  #bedfiles path
  my $bedFileDirectory = dirname($fasta_file_path);
  my @parsedHeader = parse_fasta_header($fasta_file_path);
  my $unchecked_accession = $parsedHeader[0];
  print $unchecked_accession;
  my $scientificName = $parsedHeader[1];
  print $scientificName;
  my $accession = valid_ncbi_accession($unchecked_accession);
  #check program dependencies
  my $module_path = $INC{"Bio/ViennaNGS/UCSC.pm"};
  my $template_path = $module_path;
  $template_path =~ s/UCSC.pm/template\//;
  die("ERROR [Bio::ViennaNGS::UCSC] template directory not found\n") unless (-d $template_path);
  my $faToTwoBit_path = can_run('faToTwoBit') or die 'ERROR [ViennaNGS::UCSC] faToTwoBit is not installed!';

  #create assembly hub directory structure
  my $assembly_hub_name = "assemblyHub";
  my $assembly_hub_directory = $assembly_hub_destination_path . $assembly_hub_name;
  my $genome_assembly_name = "$accession";
  my $genome_assembly_directory = $assembly_hub_directory ."/" . $genome_assembly_name;
  mkdir $assembly_hub_directory;
  mkdir $genome_assembly_directory;

  #2-bit fasta file conversion
  my $twoBitFastaFilePath = $genome_assembly_directory ."/" . "$accession" . ".2bit";
  my $full_path = can_run('faToTwoBit') or die 'faToTwoBit is not installed!';
  my $fastaToTwobit_cmd = $faToTwoBit_path . " " . $fasta_file_path . " " . $twoBitFastaFilePath;
  system($fastaToTwobit_cmd);

  #template definition
  my $template = Template->new({
                                INCLUDE_PATH => ["$template_path"],
                                RELATIVE=>1,
  });

  #construct hub.txt
  my $hubtxt_path = $assembly_hub_directory . "/hub.txt";
  my $hubtxt_file = 'hub.txt';
  my $hubtxt_vars =
    {
     hubName => "$accession",
     shortLabel => "$accession",
     longLabel => "$accession",
     genomesFile => "genome.txt",
     email => "email",
     descriptionURL => $base_URL . "description.html"
    };
  $template->process($hubtxt_file,$hubtxt_vars,$hubtxt_path) || die "Template process failed: ", $template->error(), "\n";

  #construct genome.txt
  my $genometxt_path = $assembly_hub_directory . "/genome.txt";
  my $genometxt_file = 'genome.txt';
  my $genometxt_vars =
    {
     genome => "$accession",
     trackDb => "$accession/trackDb.txt",
     groups => "$accession/groups.txt",
     description => "$accession",
     twoBitPath => "$accession/$accession.2bit",
     organism => "organism",
     defaultPos => $accession . ":1-1,000",
     orderKey => "10",
     scientificName => "$scientificName",
     htmlPath => "$accession/description.html"
    };
  $template->process($genometxt_file,$genometxt_vars,$genometxt_path) || die "Template process failed: ", $template->error(), "\n";

  #construct description.html
  my $description_html_path = $genome_assembly_directory . "/description.html";
  my $description_html_file = 'description.html';
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
  $template->process($description_html_file,$description_html_vars,$description_html_path) || die "Template process failed: ", $template->error(), "\n";

  my $groups = make_group("annotation", "Annotation", "1", "0");

  #construct group.txt
  my $group_txt_path = $genome_assembly_directory . "/groups.txt";
  my $group_txt_file = 'groups.txt';
  my $group_txt_vars = 
    {
     groups  => "$groups",
    };
  $template->process($group_txt_file,$group_txt_vars,$group_txt_path) || die "Template process failed: ", $template->error(), "\n";

  my @trackfiles = retrieve_tracks("/home/mescalin/egg/current/Projects/Perl/assemblyhub_test/", $base_URL, $assembly_hub_name);
  my $tracksList;
  foreach my $track (@trackfiles){
    my $trackString = make_track(@$track);
    $tracksList .= $trackString;
  }

  #construct trackDb.txt
  my $trackDb_txt_path = $genome_assembly_directory . "/trackDb.txt";
  my $trackDb_txt_file = 'trackDb.txt';
  my $trackDb_txt_vars = 
    {
     tracks  => "$tracksList"
    };
  $template->process($trackDb_txt_file,$trackDb_txt_vars,$trackDb_txt_path) || die "Template process failed: ", $template->error(), "\n";
}

sub retrieve_tracks{
  my ($directoryPath,$base_URL,$assembly_hub_name) = @_;
  my $currentDirectory = getcwd;
  chdir $directoryPath or die $!;
  my @trackfiles = <*.bb>;
  my @tracks;
  foreach my $trackfile (@trackfiles){
    my $filename = $trackfile;
    $filename =~ s/.bb$//;
    my @filenameSplit = split(/\./, $filename);
    my $accession = $filenameSplit[0];
    my $tag = $filenameSplit[1];
    my $id = lc($tag);
    my $track = "refseq_" . $id;
    my $bigDataUrl = $base_URL . $assembly_hub_name ."/". $trackfile;
    my $shortLabel = "RefSeq " . $tag;
    my $longLabel = "RefSeq " . $tag;
    my $type = "bigBed 12 .";
    my $autoScale = "off";
    my $bedNameLabel = "Gene Id";
    my $searchIndex = "name";
    my $colorByStrand = "100,205,255 55,155,205";
    my $visibility = "pack";
    my $group = "annotation";
    my $priority = "10";
    my @track = ($tag,$track, $bigDataUrl,$shortLabel,$longLabel,$type,$autoScale,$bedNameLabel,$searchIndex,$colorByStrand,$visibility,$group,$priority);
    my $trackreference = \@track;
    push(@tracks, $trackreference);
  }
  chdir $currentDirectory or die $!;
  return @tracks;
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
  print $fastaheader;
  close $file;
  #>gi|556503834|ref|NC_000913.3| Escherichia coli str. K-12 substr. MG1655
  my @headerfields = split(/\|/, $fastaheader);
  my $accession = $headerfields[3];
  my $scientificName = $headerfields[4];
  return ($accession,$scientificName);
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

=head2 EXPORT

Routines: 
  make_assembly_hub

Variables:
  none


=head3 make_assembly_hub()

Documentation for this routine goes here.

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
