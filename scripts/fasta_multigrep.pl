#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 17:07:55 michl>
#
# Extract individual sequences from a multi-Fasta file
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
use Cwd;
use Path::Class;
use Bio::ViennaNGS::Fasta;
use Carp;

my ($fastaO,$grep,$line,$cwd,$f,$seq,$theid, $intervals);
my @fids = ();
my $fa_in  = undef;


Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("l|list=s" => \$grep,
					   "f|fa=s"   => \$fa_in,
					   "man"      => sub{pod2usage(-verbose => 2)},
					   "help|h"   => sub{pod2usage(1)}
					   );

unless(-f $fa_in){
  warn "Could not find input FASTA input file provided via --fa option";
  pod2usage(-verbose => 0);
}
$cwd = getcwd();

unless ($fa_in =~ /^\// || $fa_in =~ /^\.\//){$fa_in = file($cwd,$fa_in);}

$fastaO = Bio::ViennaNGS::Fasta->new(fasta=>"$fa_in");

my @list = split /:/, $grep;
foreach my $id (@list){
  croak "ERROR: ID $id not found in Fasta file $fa_in"
    unless ($fastaO->has_sequid($id));
 # print ">$id\n";
  my $seq = $fastaO->fasta->get_Seq_by_id($id);

  print join "\n", (">$id", unpack "(a70)*",$seq->seq);
  print "\n";
}

__END__


=head1 NAME

fasta_multigrep.pl - Extract individual sequences from a multi Fasta file

=head1 SYNOPSIS

fasta_multigrep.pl [--list I<LIST>] [--fa I<FILE>] [options]

=head1 DESCRIPTION

This script extracts individual sequences from a multi Fasta
file. Output is written to STDOUT.

=head1 OPTIONS

=over

=item B<-l|--list>

Colon (:) separated list of Fasta IDs to extract (mandatory)

=item B<-f|--fa>

Input file in Fasta format (mandatory)

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

