# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-10 19:01:14 michl>

=head1 NAME

Bio::ViennaNGS::Fasta - Moose wrapper for Bio::DB::Fasta

=head1 SYNOPSIS

  use Bio::ViennaNGS::Fasta;

  my $f = Bio::ViennaNGS::Fasta->new(fasta => "data/foo.fa", );

  # get all FASTA IDs
  my @ids = $f->fastaids;

  # get a reference to a hash of Bio::PrimarySeq::Fasta objects whose
  # keys are the Fasta IDs seen in the input file
  my $ps = $f->primaryseqH;

  # get the strand-specific genomic sequence for a certain Fasta ID
  my $id = "chr1";
  my $start = 287;
  my $end = 1289;
  my $strand = "+";
  my $seq = $f->stranded_subsequence($id,$start,$end,$strand);

=head1 DESCRIPTION

L<Bio::ViennaNGS::Fasta> provides a L<Moose> interface to
L<Bio::DB::Fasta>, spiced up with a few convenience methods for easy
sequence data retrieval.

=head2 ATTRIBUTES

=over 3

=item fasta (required)

Upcon object construction, this attribute expects an input fasta file,
which is transparently coerced into a L<Bio::DB::Fasta> object and
hitherto available via the C<fasta> attribute.

=item fastaids (auto-computed)

Arrary reference to the Fasta IDs found in the input file

=item primaryseqH (auto-computed)

Hash reference to L<Bio::PrimarySeq::Fasta> objects whose keys are the
Fasta IDs found in the input file

=back

=cut

package Bio::ViennaNGS::Fasta;

use Bio::ViennaNGS;
use Moose;
use Bio::ViennaNGS::Subtypes;
use Bio::Perl;
use Carp;
use Data::Dumper;
use namespace::autoclean;
use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");

has 'fasta' => (
		is => 'ro',
		isa => 'Bio::ViennaNGS::MyFasta',
		required => 1,
		predicate => 'has_fasta',
		coerce => 1,
	    );


has 'fastaids' => (
		   is => 'rw',
		   isa => 'ArrayRef',
		   predicate => 'has_ids',
		   init_arg => undef,
		  );

has 'primaryseqH' => (
		     is => 'rw',
		     isa => 'HashRef',
		     predicate => 'has_primaryseq',
		     init_arg => undef,
		    );


sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  $self->fastaids([$self->fasta->ids]);
  confess "ERROR [$this_function] \$self->fastsids not available"
    unless ($self->has_ids);
  my %ps = ();
  foreach my $id (@{$self->fastaids}){
    $ps{$id} = $self->fasta->get_Seq_by_id($id);
  }
  $self->primaryseqH(\%ps);
}

=head2 METHODS

=over 2

=item stranded_subsequence

Title : stranded_subsequence

Usage : C<$obj-E<gt>stranded_subsequence($id,$start,$end,$strand)>

Function : Returns the DNA/RNA sequence for ID C<$id> from C<$start>
           to C<$end>. Internally, sequence data is retrieved from
           C<$self-E<gt>$primaryseqH> HashRef to
           L<Bio::PrimarySeqI>/L<Bio::PrimarySeq::Fasta> objects.

Args : C<$id> is the Fasta ID to retrieve sequence data from,
        C<$start> and C<$end> are (1-based) start and end coordinates
        of the requested interval, where C<$start> must be <= C<$end>,
        and C<$strand> is 1 or -1 for [+] or [-] strand, respectively.

Returns : A string.

=item has_sequid

Title : has_sequid 

Usage : C<$obj-E<gt>has_seqid($id)>

Function : Checks whether the current object contains Fasta ID C<$id>.

Args : C<$id> if the Fasta ID to check for.

Returns : 1 if ID C<$id> was found, 0 else.

=back

=cut

sub stranded_subsequence {
  my ($self,$id,$start,$end,$strand) = @_;
  my ($this_function,$seq,$rc,$p,$obj);
  $this_function = (caller(0))[3];
  confess "ERROR [$this_function] start corrdinate must be <= end coordinate"
    unless ($start <= $end);
  confess "ERROR [$this_function] Id $id not found in input Fasta file"
    unless ($self->has_sequid($id));
  $seq = ${$self->primaryseqH}{$id}->subseq($start => $end);
  if ($strand eq '-1' || $strand eq '-') {
    $rc = revcom($seq);
    $seq = $rc->seq();
  }
  #print "id:$id\nstart:$start\nend:$end\n";
  return $seq;
}


sub has_sequid {
  my ($self,$id) = @_;
  return exists ${$self->primaryseqH}{$id} ? 1 : 0;
}


=head1 DEPENDENCIES

=over

=item L<Bio::Perl> >= 1.00690001

=item L<Bio::DB::Fasta>

=item L<Moose>

=item L<Carp>

=item L<namespace::autoclean>

=back

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::DB::Fasta>

=back

=head1 AUTHOR

Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2017 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

1;
