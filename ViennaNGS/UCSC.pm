# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-22 17:01:47 mtw>

package Bio::ViennaNGS::UCSC;

use Exporter;
use version; our $VERSION = qv('0.01');
use strict;
use warnings;
use Template;
use File::Which qw(which where);

our @ISA = qw(Exporter);

our @EXPORT_OK = qw( make_assembly_hub  );

our @EXPORT = ();

sub make_assembly_hub{
  my ($fasta_file_path,$assembly_hub_destination_path,$base_URL,$log_path) = @_;
  #check arguments
  die("ERROR [ViennaNGS::UCSC] $fasta_file_path does not exist\n") unless (-e $fasta_file_path);
  die("ERROR [ViennaNGS::UCSC] $assembly_hub_destination_path does not exist\n") unless (-d $assembly_hub_destination_path);
  die("ERROR [ViennaNGS::UCSC] $base_URL is not set\n") if ($base_URL == "");
  die("ERROR [ViennaNGS::UCSC] $log_path does not exist\n") unless (-e $log_path);

  #check program dependencies
  my $faToTwoBit_path = which('faToTwoBit');
  die("ERROR [ViennaNGS::UCSC] faToTwoBit does not exist in $PATH\n") unless (-e $faToTwoBit_path);
  my $gff2bed_path = which('gff2bed');
  die("ERROR [ViennaNGS::UCSC] gff2bed does not exist in $PATH\n") unless (-e $gff2bed);

  #create assembly hub directory structure
  my $assembly_hub_name = "Test";
  my $assembly_hub_directory = $assembly_hub_destination_path . $assembly_hub_name;
  my $genome_assembly_name = "genomeAssembly";
  my $genome_assembly_directory = $assembly_hub_directory ."/" . $genome_assembly_name;
  mkdir $assembly_hub_directory;
  mkdir $genome_assembly_directory;

  #2-bit fasta file conversion
  my $2bit_fasta_file_path = $genome_assembly_directory ."/" . "genomeAssembly" . ".2bit"
  my $fastaToTwobit_cmd =  $faToTwoBit_path . " " . $fasta_file_path . " " . $2bit_fasta_file_path;
  print STDERR ">> " . $fastaToTwobit_cmd . "\n";
  system($fastaToTwobit_cmd);
}

1;
__END__

=head1 NAME

ViennaNGS::UCSC - Perl extension for easy UCSC Genome Browser
integration.

=head1 SYNOPSIS

  use ViennaNGS::UCSC;

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
Florian Eggenhofer, E<lt>egg@tbi.univie.ac.atE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


=cut
