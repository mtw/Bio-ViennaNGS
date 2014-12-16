# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-17 00:11:48 mtw>

package Bio::ViennaNGS::BamStatSingle;

use 5.12.0;
use version; our $VERSION = qv('0.12_05');
use Bio::DB::Sam 1.39;
use Moose;
use Carp;
use Data::Dumper;
use namespace::autoclean;

has 'bam' => (
	      is => 'rw',
	      isa => 'Str',
	      required => 1,
	      predicate => 'has_bam',
	      trigger => \&_parse_bam,
	     );

has 'data' => (
	       is => 'ro',
	       isa => 'HashRef',
	       predicate => 'has_data',	       
	      );

sub _parse_bam {
  my ($self) = @_;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] Attribute 'bam' not found $!"
    unless ($self->has_bam);
  print "parsing BAM file $self->bam ...\n";

}

__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::BamStatSingle - Moose interface to BAM mapping statistics

=head1 SYNOPSIS

  use Bio::ViennaNGS::BamStatSingle;

  my $bss1 = Bio::ViennaNGS::BamStatSingle->new(bam => "path/to/file.bam");

=head1 DESCRIPTION

This module provides a L<Moose> interface to the mapping statistics of
a single BAM file. It builds on L<Bio::DB::Sam> and 

=head1 SEE ALSO

=over 

=item L<Bio::ViennaNGS>

=back

=head1 AUTHORS

=over 

=item Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=item Michael T. Wolfinger  E<lt>michael@wolfinger.euE<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
