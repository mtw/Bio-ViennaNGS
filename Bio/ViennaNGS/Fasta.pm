package Bio::ViennaNGS::Fasta;

use 5.12.0;
use version; our $VERSION = qv('0.02_01');
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

__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::Fasta - Moose warapper for Bio::DB::Fasta

=head1 SYNOPSIS

  use Bio::ViennaNGS::Fasta;

  my $f = Bio::ViennaNGS::Fasta->new( fa => "data/foo.fa", );

  # get all FASTA IDs
  my @ids = $f->fastaids;

=head1 DESCRIPTION

This module provides a L<Moose> interface to L<Bio::DB::Fasta>.

=head1 EXPORT

None by default.


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
