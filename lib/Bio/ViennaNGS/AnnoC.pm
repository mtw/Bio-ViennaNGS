# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:30:51 mtw>

package Bio::ViennaNGS::AnnoC;

use 5.12.0;
use version; our $VERSION = qv('0.12_07');
use Bio::ViennaNGS qw(sortbed);
use Bio::Tools::GFF;
use Path::Class;
use Carp;
use Moose;
use namespace::autoclean;

has 'accession' => (
		    is => 'rw',
		    isa => 'Str',
		    predicate => 'has_accession',
		   );

has 'features' => (
		   is => 'ro',
		   isa => 'HashRef',
		   predicate => 'has_features',
		   default => sub { {} },
		  );

has 'nr_features' => (
		     is => 'ro',
		     isa => 'Int',
		     builder => '_get_nr_of_features',
		     lazy => 1,
		     );

has 'featstat' => (
		   is => 'ro',
		   isa => 'HashRef',
		   builder => '_set_featstat',
		   predicate => 'has_featstat',
		   lazy => 1,
		  );

before 'featstat' => sub {
  my $self = shift;
  $self->_get_nr_of_features();
};

sub _set_featstat {
  my $self = shift;
  my $this_function = (caller(0))[3];
  my %fs = ();
  confess "ERROR [$this_function] \$self->features not available"
    unless ($self->has_features);
  $fs{total} = 0;
  $fs{origin} = "$this_function ".$VERSION;
  $fs{count} = $self->nr_features;
  foreach my $uid ( keys %{$self->features} ){
    my $gbkey = ${$self->features}{$uid}->{gbkey};
    $fs{total} += 1;
    unless (exists $fs{$gbkey}){
      $fs{$gbkey} = 0;
    }
    $fs{$gbkey} += 1;
  }
  return \%fs;
}

sub _get_nr_of_features {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->features not available"
    unless ($self->has_features);
  return (keys %{$self->features});
}

sub parse_gff {
  my ($self,$in_file) = @_;
  my ($i,$gffio,$header,$f,$gbkey);
  my $this_function = (caller(0))[3];

  $gffio = Bio::Tools::GFF->new(-file         => $in_file,
				-gff_version  => 3,
			       );
  $gffio->ignore_sequence(1);
  if ($header = $gffio->next_segment() ){
    $self->accession( $header->display_id() );
  }
  else{ carp "ERROR [$this_function] Cannot parse GFF header\n" }

  while($f = $gffio->next_feature()) {
    my ($uid,$feat_name);
    my @name = my @id = my @gbkeys = ();

    next if ($f->primary_tag() eq "exon");

    # 1) determine gbkey of the current feature
    @gbkeys = $f->get_tag_values("gbkey");
    $gbkey  = $gbkeys[0];

    # 2) get a unique ID for each feature
    if ($f->has_tag('ID')){
      @id = $f->get_tag_values('ID');
      $uid = $id[0]; # ID=id101
    }
    else {
      croak "ERROR [$this_function] Feature '$gbkey' at pos.\
             $f->start does not have \'ID\' attribute\n";
    }

    # 3) assign parent's unique ID in case a parent record exists
    if ($f->has_tag('Parent')){
      @id = $f->get_tag_values('Parent');
      $uid = $id[0]; # ID=id101
    }

    # 4) find a name for the current feature, use 'Name' or 'ID' attribute
    if ($f->has_tag('Name')){
      @name = $f->get_tag_values('Name');
      $feat_name = $name[0];
    }
    elsif ($f->has_tag('ID')){
      @id = $f->get_tag_values('ID');
      $feat_name = $id[0]; # ID=id101, use ID as feature name
    }
    else {
      croak "ERROR [$this_function] Cannot set name for feature \
              $f->gbkey at pos. $f->start\n";
    }

    unless (exists ${$self->features}{$uid}) { # gene / ribosome_entry_site / etc.
      ${$self->features}{$uid}->{start}   = $f->start;
      ${$self->features}{$uid}->{end}     = $f->end;
      ${$self->features}{$uid}->{strand}  = $f->strand;
      ${$self->features}{$uid}->{length}  = $f->length;
      ${$self->features}{$uid}->{seqid}   = $f->seq_id;
      ${$self->features}{$uid}->{score}   = $f->score || 0;
      ${$self->features}{$uid}->{gbkey}   = $gbkey;
      ${$self->features}{$uid}->{name}    = $feat_name;
      ${$self->features}{$uid}->{uid}     = $uid;
    }
    else { # CDS / tRNA / rRNA / etc.
      ${$self->features}{$uid}->{gbkey} = $gbkey;  # gbkey for tRNA/ rRNA/ CDS etc
    }
  }
  $gffio->close();
}

sub features2bed {
 my ($self,$gbkey,$dest,$bn,$log) = @_;
 my ($chrom,$chrom_start,$chrom_end,$name,$score,$strand,$thick_start);
 my ($thick_end,$reserved,$block_count,$block_sizes,$block_starts);
 my @ft = ();
 my $this_function = (caller(0))[3];

 croak "ERROR [$this_function] $self->features not available"
   unless ($self->has_features);
 croak "ERROT [$this_function] $self->featstat not available"
   unless ($self->has_featstat);
 croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);
 if (defined $log){open(LOG, ">>", $log) or croak $!;}

 if (defined $gbkey){ # dump features of just one genbank key
   confess "ERROR [$this_function] genbank key \'$gbkey\' N/A in hash "
     unless (exists ${$self->featstat}{$gbkey});
   $ft[0] = $gbkey;
 }
 else{ # dump features for all genbank keys
   foreach my $gbk (keys %{$self->featstat}) {
     next if ($gbk eq 'total' || $gbk eq 'Src' || $gbk eq 'accession' ||
	      $gbk eq 'origin' || $gbk eq 'count');
     push @ft,$gbk;
   }
 }

 foreach my $f (@ft){
   my $bedname   = file($dest,"$bn.$f.bed");
   my $bedname_u = file($dest,"$bn.$f.u.bed");
   open (BEDOUT, "> $bedname_u") or croak $!;

   # dump unsorted gene annotation from DS to BED12
   foreach my $uid (keys %{$self->features}){
      next unless (${$self->features}{$uid}->{gbkey} eq $f);
      my @bedline = ();
      $chrom        = ${$self->features}{$uid}->{seqid};
      $chrom_start  = ${$self->features}{$uid}->{start};
      $chrom_start--; # BED is 0-based
      $chrom_end    = ${$self->features}{$uid}->{end};
      $name         = ${$self->features}{$uid}->{name};
      $score        = ${$self->features}{$uid}->{score};
      $strand       = ${$self->features}{$uid}->{strand} == -1 ? '-' : '+'; #default to +
      $thick_start  = $chrom_start;
      $thick_end    = $chrom_end;
      $reserved     = 0; 
      $block_count  = 1;
      $block_sizes  = ${$self->features}{$uid}->{length}.",";
      $block_starts = "0,";
      @bedline = join ("\t", ($chrom,$chrom_start,$chrom_end,
			      $name,$score,$strand,$thick_start,
			      $thick_end,$reserved,$block_count,
			      $block_sizes, $block_starts));
      print BEDOUT "@bedline\n";
    }
   close (BEDOUT);

   sortbed($bedname_u,".",$bedname,1,undef);  # sort bed file

 } # end foreach
 if (defined $log){close(LOG)};
}

sub feature_summary {
  my ($self, $dest) = @_;
  my ($fn,$fh);
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] $dest does not exist\n"
    unless (-d $dest);
  croak "ERROR [$this_function] $self->accession not available\n"
    unless ($self->has_accession);

  $fn = dir($dest,$self->accession.".summary.txt");
  open $fh, ">", $fn or croak $!;

  print $fh "Accession\t ".$self->accession."\n";
  print $fh "Origin   \t ${$self->featstat}{origin}\n";
  foreach my $ft (sort keys %{$self->featstat}){
    next if ($ft =~ /total/ || $ft =~ /accession/ || $ft =~ /origin/);
    print $fh "$ft\t${$self->featstat}{$ft}\n";
  }
  print $fh "Total\t${$self->featstat}{total}\n";
  close $fh;
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 NAME

Bio::ViennaNGS::AnnoC - Object-oriented interface for storing and
converting biological sequence annotation formats

=head1 SYNOPSIS

  use Bio::ViennaNGS::AnnoC;

  my $obj = Bio::ViennaNGS::AnnoC->new();

  # parse GFF3 file to internal data straucture
  $obj->parse_gff($gff3_file);

  # compute summary of parsed annotation
  $obj->featstat;

  # dump feature summary to file
  $obj->feature_summary($dest);

  # dump all tRNAs contained in data structure as BED12
  $obj->features2bed("tRNA",$dest,$bn,$log)

=head1 DESCRIPTION

This module provides an object-oriented interface for storing and
converting biological sequence annotation data. Based on the C<Moose>
object system, it maintains a central data structure which is curently
designed to represent simple, non-spliced (ie single-exon) annotation
data. Future versions of the module will account for more generic
scenarios, including spliced isoforms.

=head1 METHODS

=over

=item parse_gff

 Title   : parse_gff
 Usage   : $obj->parse_gff($gff3_file);
 Function: Parses GFF3 annotation files of non-spliced genomes into
           C<$self->features>
 Args    : The full path to a GFF3 file
 Returns :
 Notes   : The GFF3 specification is available at
           L<http://www.sequenceontology.org/resources/gff3.html>.
           This routine has been tested with NCBI bacteria GFF3
           annotation.

=item feature_summary

 Title   : feature_summary
 Usage   : $obj->feature_summary($dest);
 Function: Generate a summary file for all features present in
           C<$self->features> 
 Args    : Full output path for summary.txt file
 Returns :


=item features2bed

 Title   : features2bed
 Usage   : $obj->features2bed($feature,$workdir,$bn,$log);
 Function: Dumps genomic features from C<$self->features> hash to a
           BED12 file.
 Args    : C<$gbkey> can be either a string corresponding to a
           genbank key in C<$self->featstat> or C<undef>. If defined,
           only features of the speficied key will be dumped to a single
           BED12 file. If C<$gbkey> is C<undef>, BED12 files will be
           generated for each type present in C<$self->featstat>.
           C<$dest> is the output directory and C<$bn> the basename for
           all output files. C<$log> is either be the full path to a
           logfile or C<undef>.
 Returns  :

=back

=head1 DEPENDENCIES

=over 4

=item L<Bio::Tools::GFF>

=item L<IPC::Cmd>

=item L<Path::Class>

=item L<Carp>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.


=cut
