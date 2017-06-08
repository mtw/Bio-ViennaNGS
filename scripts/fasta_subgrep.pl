#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 12:38:25 michl>
#
# Extract subsequences from a (multi) Fasta file
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2017 Michael T. Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# *  This program is free software: you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *  This copyright notice MUST APPEAR in all copies of the script!
# ***********************************************************************

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Data::Dumper;
use Path::Class;
use Bio::ViennaNGS::Fasta;
use Carp;
use Cwd;
use Data::Dumper;

my ($fastaO,$cwd,$seq);
my $fa_in = undef;
my $start = undef;
my $end = undef;
my $id = undef;
my $strand = "+";

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions('f|fa=s'     => \$fa_in,
					   's|start=s'  => \$start,
					   'e|end=s'    => \$end,
					   'id=s'       => \$id,
					   'strand=s'   => \$strand);

unless (-f $fa_in){
  warn "Could not find input FASTA input file provided via -f|--fa option";
  pod2usage(-verbose => 0);
}
$cwd = getcwd();
unless ($fa_in  =~ /^\// || $fa_in  =~ /^\.\//){$fa_in = file($cwd,$fa_in);}
croak "ERROR: Option -s|--start must be provided"
  unless (defined($start));
croak "ERROR: Option -e|--end must be provided"
  unless (defined($end));
$strand eq '-' ? ($strand = '-') : ($strand = '+') ;

$fastaO = Bio::ViennaNGS::Fasta->new(fasta=>"$fa_in");
#print Dumper($fastaO);
if(scalar(@{$fastaO->fastaids}) > 1){ # multi-Fasta file
  croak "ERROR: Option --id must be provided for multi-Fasta input files"
    unless (defined($id));
}
else { # single-Fasta file
  carp "WARNING: The ID provided ($id) does not match the single ID presentin the input Fasta file (${$fastaO->fastaids}[0]). Using ${$fastaO->fastaids}[0] !"
    unless ($id eq ${$fastaO->fastaids}[0]);
    $id = ${$fastaO->fastaids}[0];
}

confess "[ERROR]: ID $id not found in Fasta file $fa_in"
    unless ($fastaO->has_sequid($id));
$seq = $fastaO->stranded_subsequence($id,
				     $start,
				     $end,
				     $strand);
print ">$id\n";
print  join "\n", (unpack "(a70)*",$seq);
print "\n";

__END__


=head1 NAME

fasta_subgrep.pl - Extract subsequence from a (multi) Fasta file

=head1 SYNOPSIS

fasta_subgrep.pl [--fa I<FILE>] [-s I<INT>] [-e I<INT>] [-id I<STRING>] [--strand I<+/->]

=head1 DESCRIPTION

Extract a subsequence from a (multi) Fasta file. Interval coordinates
are passed as start and end coordinates via the C<-s> and C<-e>
options, respectively. Output is written to STDOUT.

=head1 OPTIONS

=over

=item B<-f|--fa>

Input file in Fasta format (mandatory)

=item B<-s|--start>

Start of the sequence interval to extract

=item B<-e|--end>

End of the sequence interval to extract

=item B<--id>

Fasta ID of the sequence to extract from. This is only required for
multi Fasta input files.

=item B<--strand>

Specify the strand to extract sequence data from. Allowed arguments
are C<+> amd C<->. In the latter case, retrieve reverse complement of
the sequence interval from C<start> to C<end>.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

