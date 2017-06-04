#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-03-15 17:27:44 michl>
#
# Extract sequences from a multi-Fasta file
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
use Bio::DB::Fasta;

#use Bio::SeqUtils;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($fastaO,$grep,$line,$cwd,$f,$seq,$theid, $intervals);
my @fids = ();
my $fa_in  = undef;


Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("list=s"   => \$grep,
					   "fa=s"     => \$fa_in,
					   "man"        => sub{pod2usage(-verbose => 2)},
					   "help|h"     => sub{pod2usage(1)}
					   );

unless(-f $fa_in){
  warn "Could not find input FASTA input file provided via ---fa option";
  pod2usage(-verbose => 0);
}
$cwd = getcwd();

unless ($fa_in =~ /^\// || $fa_in =~ /^\.\//){$fa_in = file($cwd,$fa_in);}

$fastaO = Bio::DB::Fasta->new($fa_in);
#@fids = $fastaO->fastaids; # get all FASTA IDs
#print Dumper(@fids);
#print Dumper($fastaO);

my @list = split /:/, $grep;
foreach my $id (@list){
  print ">$id\n";
  my $seq = $fastaO->get_Seq_by_id($id);
  unless(defined $seq){
    die "ID $id not found in file $fa_in";
  }
  print join "\n", (unpack "(a70)*",$seq->seq);
  print "\n";
}

__END__


=head1 NAME

grep_multifa.pl - Extract individual sequences from a multi Fasta file

=head1 SYNOPSIS

grep_multifa.pl [--list I<LIST>] [--fa I<FILE>] [--aa] [options]

=head1 DESCRIPTION

This script extracts individual sequences from a multi Fasta file

=head1 OPTIONS

=over

=item B<--list>

Colon (:) separated list of Fasta IDs to extract

=item B<--fa>

Input file in Fasta format (mandatory)

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

