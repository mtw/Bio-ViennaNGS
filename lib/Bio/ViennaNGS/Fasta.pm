package Bio::ViennaNGS::Fasta;

use 5.12.0;
use version; our $VERSION = qv('0.11');
use Bio::Perl 1.00690001;
use Bio::DB::Fasta;
use Moose;
use Carp;
use Data::Dumper;
use namespace::autoclean;

has 'fa' => (
	     is => 'rw',
	     isa => 'Str',
	     required => 1,
	     predicate => 'has_fa',
	    );

has 'fastadb' => (
		  is => 'rw',
		  isa => 'Bio::DB::Fasta',
		  builder => '_get_fastadb',
		  predicate => 'has_db',
		  lazy => 1,
		 );

has 'fastaids' => (
		   is => 'ro',
		   isa => 'ArrayRef',
		   builder => '_get_fastaids',
		   predicate => 'has_ids',
		   lazy => 1,
		  );

has 'primaryseq' => (
		     is => 'ro',
		     isa => 'HashRef',
		     builder => '_get_primaryseq',
		     lazy => 1,
		    );

before 'primaryseq' => sub{
  my $self = shift;
  $self->fastaids;
};

before 'fastaids' => sub{
  my $self = shift;
  $self->fastadb;
};

sub _get_fastadb {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] Fasta input not available"
    unless (-f $self->fa);
  my $db =  Bio::DB::Fasta->new($self->fa) or croak $!;
  return $db;
}

sub _get_fastaids {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] Attribute 'fastadb' not found $!"
    unless ($self->has_db);
  my $db = $self->fastadb or croak $!;
  my @ids = $db->ids or croak $!;
  return \@ids;
}

sub _get_primaryseq {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] Attribute 'fastaids' not found $!"
    unless ($self->has_ids);
  my %fobj = ();
  my $db = $self->fastadb or croak $!;
  foreach my $id (@{$self->fastaids}) {
    $fobj{$id} = $db->get_Seq_by_id($id); # Bio::PrimarySeq::Fasta object
  }
  return \%fobj;
}

# stranded_subsequence ($id,$start,$stop,$strand)
# retrieve RNA/DNA sequence from a Bio::PrimarySeqI /
# Bio::PrimarySeq::Fasta object
sub stranded_subsequence {
  my ($self,$id,$start,$end,$strand) = @_;
  my ($this_function,$seq,$rc,$p,$obj);
  $this_function = (caller(0))[3];
  my @dummy = $self->fastaids;
  confess "ERROR [$this_function] Attribute 'fastaids' not found $!"
    unless ($self->has_ids);
  $p = $self->primaryseq; # Hash of Bio::PrimarySeq::Fasta objects
  confess "ERROR [$this_function] Fasta ID $id not found in hash $!"
    unless (exists $$p{$id});
  $obj = $$p{$id};
  $seq = $obj->subseq($start => $end);
  if ($strand eq '-1' || $strand eq '-') {
    $rc = revcom($seq);
    $seq = $rc->seq();
  }
  # print "id:$id\nstart:$start\nend:$end\n";
  return $seq;
}

sub has_sequid {
  my ($self,$id) = @_;
  my $ids = $self->fastaids;
  #$i = grep{$_ eq $id}@{$ids} ? 1 : 0;
  for my $j (@$ids){
    if ($id eq $j){return 1;}
    else {return 0;}
  }
  return -1;
}

__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::Fasta - Moose wrapper for Bio::DB::Fasta

=head1 SYNOPSIS

  use Bio::ViennaNGS::Fasta;

  my $f = Bio::ViennaNGS::Fasta->new( fa => "data/foo.fa", );

  # get all FASTA IDs
  my @ids = $f->fastaids;

  # get a reference to a hash of Bio::PrimarySeq::Fasta objects whose
  # keys are the Fasta IDs
  my $ps = $f->primaryseq;

  # get the strand-specific genomic sequence for a certain Fasta ID
  my $id = "chr1";
  my $start = 287;
  my $end = 1289;
  my $strand = "+";
  my $seq = $foo->stranded_subsequence($id,$start,$end,$strand);

=head1 DESCRIPTION

This module provides a L<Moose> interface to L<Bio::DB::Fasta>.

=head1 METHODS

=over

=item stranded_subsequence

  Title   : stranded_subsequence
  Usage   : $obj->stranded_subsequence($id,$start,$end,$strand)
  Function: Returns the DNA/RNA sequence for C<$id> from
            C<$start> to C<$end>.
  Args    : C<$id> is the Fasta ID (a L<Bio::PrimarySeq::Fasta> object). 
            C<$start> and C<$end>  should be self-explnatory, C<$strand>
            is 1 or -1 for [+] or [-] strand, respectively
  Returns : A string

=back

=head1 DEPENDENCIES

=over 5

=item L<Bio::Perl> >= 1.00690001

=item L<Bio::DB::Fasta>

=item L<Moose>

=item L<Carp>

=item L<namespace::autoclean>

=back

=head1 SEE ALSO

=over 2

=item L<Bio::ViennaNGS>

=item L<Bio::DB::Fasta>

=back

=head1 AUTHOR

Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
