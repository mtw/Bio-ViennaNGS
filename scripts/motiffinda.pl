#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-09 23:58:37 mtw>
#
# Find motifs in annotated sequences. The motif is provided as regex
# via the command line
#
# usage: motiffinda.pl -motif <REXGEX> -gff <GFFFILE> -fa <FASTAFILE>
#                      -offset <OFFSET> -inframe
#
# The motif must be privided within braces, because we use $1 to obtain
# the position of the motif within the query sequence.
#
# ./motiffinda.pl -motif "(ACATG\w{4,13}ACATG)" -gff NC_000913.gff
#                 -fa NC_000913.fna -offset 1
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael T. Wolfinger <michael@wolfinger.eu>
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
#

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Bio::ViennaNGS::Fasta;
use Bio::ViennaNGS::AnnoC qw(&parse_gff $feat);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($gff_in,$fa_in,$fastaO);
my $gbkey      = "CDS";
my $motif      = undef;
my $inframe    = 0;
my $offset     = 0;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("gff|g=s"    => \$gff_in,
					   "fa|f=s"     => \$fa_in,
					   "motif|m=s"  => \$motif,
					   "gbkey=s"    => \$gbkey,
					   "inframe|i"  => sub{$inframe=1},
					   "offset|o=i" => \$offset,
					   "man"        => sub{pod2usage(-verbose => 2)},
					   "help|h"     => sub{pod2usage(1)}
					   );
unless ($gff_in =~ /^\//) {$gff_in = "./".$gff_in;}
unless (-f $gff_in){
  warn "Could not find input file $gff_in given via --gff|-g option";
  pod2usage(-verbose => 0);
}

unless ($fa_in =~ /^\//) {$fa_in = "./".$fa_in;}
unless (-f $fa_in){
  warn "Could not find input file $fa_in given via --fa|-f option";
  pod2usage(-verbose => 0);
}

unless (defined $motif) {
  warn "Please provide motif as regular expression";
  pod2usage(-verbose => 0);
}

$fastaO = Bio::ViennaNGS::Fasta->new(fa=>$fa_in);

parse_gff($gff_in);
find_motifs($fastaO);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# loop over all features and check for the motif
sub find_motifs {
  my ($fo) = @_;
  my ($id,$start,$stop,$strand,$name,$seq,$length,$chr);

  foreach $id (keys %$feat){
    my $have_motif = 0;
    my @fm = ();  # AoH for storing found motifs and their position

    # TODO let user choose where to look for the mtoif 
    next unless ($$feat{$id}->{gbkey} eq $gbkey);
    $start  = $$feat{$id}->{start};
    $stop   = $$feat{$id}->{end};
    $strand = $$feat{$id}->{strand};
    $name   = $$feat{$id}->{name};
    $length = $$feat{$id}->{length};
    $chr    = $$feat{$id}->{seqid};

    # check if we can provide sequence information for this $chr
    unless ($fo->has_sequid($chr)){
      my $_ids= $fo->fastaids;
      print STDERR "\'$chr\' is not a valid Fasta ID. ";
      print STDERR "The provided Fasta file has the following IDs:\n";
      print STDERR join " | ", @$_ids;
      print STDERR "\n";
      die;
    }

    $seq    = $fo->stranded_subsequence($chr,$start,$stop,$strand);

    #  $seq    =~ s/T/U/g; # make RNA from DNA
    while ($seq =~ m/$motif/g) {
      my $p = pos($seq)-length($1)+1;
      $have_motif++;
      push (@fm, {motif=>$1,
		  pos=>$p,
		  name=>$name,
		  strand=>$strand,
		  start=>$start,
		  stop=>$stop},
	   );
    }

    for (my $i=0;$i<$have_motif;$i++) {
      next unless (($inframe == 1) && (($fm[$i]->{pos}+$offset)%3==0));
      # ATG in frame
      print ">$fm[$i]->{name}|$fm[$i]->{strand}|rel $fm[$i]->{pos}|";
      if ($fm[$i]->{strand} == 1){
	print "abs ".eval($fm[$i]->{pos}+$fm[$i]->{start}-1)."\n";
      }
      else {
	print "abs ".eval($fm[$i]->{stop}-$fm[$i]->{pos}+1)."\n";
      }
      print "$fm[$i]->{motif}\n";
    } # end for

  } # end foreach
}


__END__

=head1 NAME

motiffinda.pl - Find sequence motifs in annotated features.

=head1 SYNOPSIS

motiffinda.pl [--motif I<REGEX>] [--gff I<FILE>] [--fa I<FILE>]
[options]

=head1 DESCRIPTION

This script extracts sequence motifs from gene annotation data. The
motif can be provided as regualr expression. Optionally, only motifs
I<in frame> with the annoation are reported.

The tool returns B<all motifs> matching the search criteria as a
(mutli-)Fasta file to STDOUT. This means if teh motif is found more
than once within an annotated feature, all matches will be reported.

=head1 OPTIONS

=over

=item B<motif|m>

The motif to search for as regular expression. For technical reasons,
the regular expression must be enclosed in brackets ().

=item B<gff|g>

Genome annotation in GFF3 format

=item B<fa|f>

Reference genome in Fasta format

=item B<gbkey>

Motifs are only searched for in this feature type, aka Genbank key
(e.g. CDS or gene). Currently only one Genbank key provided via this
option will be procesed.

=item B<offset|o>

Offset for determination of frame

=item <inframe|i>

Only report motifs in current ORF


=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
