#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2013-09-16 14:03:48 mtw>
#
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael Thomas Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# * This library is free software; you can redistribute it and/or modify
# * it under the same terms as Perl itself, either Perl version 5.12.4 or,
# * at your option, any later version of Perl 5 you may have available.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# *
# ***********************************************************************

package ViennaNGS;

use Exporter;
use strict;
use warnings;
use Bio::Perl;

our @ISA = qw(Exporter);
our $VERSION = '0.01';
our @EXPORT = qw(get_stranded_subsequence);

our @EXPORT_OK = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# get_stranded_subsequence ($obj,$start,$stop,$path)
# retrieve RNA/DNA sequence from a Bio::PrimarySeqI /
# Bio::PrimarySeq::Fasta object
sub get_stranded_subsequence {
  my ($o,$start,$end,$strand) = @_;
  my $seq = $o->subseq($start => $end);
  if ($strand == -1) {
    my $rc = revcom($seq);
    $seq = $rc->seq();
  }
  #print "id:$id\nstart:$start\nend:$end\n";
  return $seq;
}


1;
__END__


=head1 NAME

ViennaNGS - Perl extension for analysis of Next-Generation Sequencing
(NGS) data.

=head1 SYNOPSIS

  use ViennaNGS;

=head1 DESCRIPTION

ViennaNGS is a collection of subroutines often used for NGS data analysis.

=head1 EXPORT

Routines: get_stranded_subsequence($obj,$start,$stop,$path)

Variables: none

=head2 get_stranded_subsequence($object,$start,$stop,$strand)

Returns the actual DNA/RNA sequence from $start to $stop. $object is a
Bio::PrimarySeq::Fasta object, which obeys the Bio::PrimarySeqI
conventions. To recover the entire raw DNA or protein sequence,
call $object->seq(). $strand is 1 or -1.

=head1 SEE ALSO

perldoc ViennaNGS::AnnoC

=head1 AUTHOR

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Michael Thomas Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.


=cut
