# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-16 00:08:43 mtw>

package Bio::ViennaNGS::AnnoC;

use version; our $VERSION = qv('0.07_01');
use Bio::Tools::GFF;
use Bio::DB::Fasta;
use IPC::Cmd qw(can_run run);
use Path::Class;
use Carp;
use Moose;
use Data::Dumper;


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
  $fs{accession} = "n/a";
  $fs{origin} = "$this_function ".$VERSION;
  $fs{count} = $self->nr_features;
  foreach my $key ( keys %{$self->features} ){
    #print Dumper (\$key);
    $fs{total} += 1;
    next if ($key eq 'accession');
    next if ($key eq 'origin');
    unless (exists $fs{$key}){
      $fs{$key} = 0;
    }
    $fs{$key} += 1;
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
  #if ($header = $gffio->next_segment() ){
  #  ${$self->featstat}{accession}= $header->display_id();
  #}
  #else{ carp "[$this_function]: Cannot parse GFF header\n" }

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
      ${$self->features}{$uid}->{start}     = $f->start;
      ${$self->features}{$uid}->{end}       = $f->end;
      ${$self->features}{$uid}->{strand}    = $f->strand;
      ${$self->features}{$uid}->{length}    = $f->length;
      ${$self->features}{$uid}->{seqid}     = $f->seq_id;
      ${$self->features}{$uid}->{score}     = $f->score || 0;
      ${$self->features}{$uid}->{gbkey}     = $gbkey;
      ${$self->features}{$uid}->{name}      = $feat_name;
      ${$self->features}{$uid}->{uid}       = $uid;
    }
    else { # CDS / tRNA / rRNA / etcx
      ${$self->features}{$uid}->{gbkey} = $gbkey;  # gbkey for tRNA/ rRNA/ CDS etc
    }
  }

  # finally generate some statistics on features present in this annotation
  #$self->_set_featstat;

  $gffio->close();
}

no Moose;

#our @ISA       = qw(Exporter);
#our @EXPORT_OK = qw(&parse_gff &feature_summary &get_fasta_ids &features2bed
#		    $feat $fstat $fastadb
#		    @fastaids
#		    %features %featsta);
#our @EXPORT    = ();
#our ($fastadb);
#our %features  = ();
our %featstat  = ();
#our $feat      = \%features;
our $fstat     = \%featstat;
#our @fastaids  = ();


# features2bed($featR,$fstatR,$feature,$dest,$bn,$log)
# Convert genome annotation features to BED12
# Returns a BED12 for each feature type in %features hash
sub features2bed {
 my ($featR,$fstatR,$feature,$dest,$bn,$log) = @_;
 my ($chrom,$chrom_start,$chrom_end,$name,$score,$strand,$thick_start);
 my ($thick_end,$reserved,$block_count,$block_sizes,$block_starts);
 my @ft = ();
 my $this_function = (caller(0))[3];
 my $bedtools = can_run('bedtools') or
   croak "ERROR [$this_function] Cannot find 'bedtools' utility";

 croak "ERROR [$this_function] undef \%features hash"
   unless (defined $featR);
 croak "ERROT [$this_function] undef \%fstat hash"
   unless (defined $fstatR);
 croak "ERROR [$this_function] $dest does not exist"
    unless (-d $dest);
  if (defined $log){open(LOG, ">>", $log) or croak $!;}

 if (defined $feature){ # dump one feature type
   confess "ERROR [$this_function] feature type \'$feature\' N/A in hash "
     unless (exists $$fstat{$feature});
   $ft[0] = $feature;
 }
 else{ # dump all feature types
   # get all feature types from in %$featR
   foreach my $gbkey (keys %$fstatR) {
     next if ($gbkey eq 'total' || $gbkey eq 'Src' ||
	      $gbkey eq 'accession' || $gbkey eq 'origin');
     push @ft,$gbkey;
   }
 }

 foreach my $f (@ft){
   my $bedname   = file($dest,"$bn.$f.bed");
   my $bedname_u = file($dest,"$bn.$f.u.bed");
   open (BEDOUT, "> $bedname_u") or croak $!;

   # dump unsorted gene annotation from DS to BED12
   foreach my $uid (keys %$featR){
      next unless ($$featR{$uid}->{gbkey} eq $f);
      my @bedline = ();
      $chrom        = $$featR{$uid}->{seqid};
      $chrom_start  = $$featR{$uid}->{start};
      $chrom_start--; # BED is 0-based
      $chrom_end    = $$featR{$uid}->{end};
      $name         = $$featR{$uid}->{name};
      $score        = $$featR{$uid}->{score};
      $strand       = $$featR{$uid}->{strand} == -1 ? '-' : '+'; #default to +
      $thick_start  = $chrom_start;
      $thick_end    = $chrom_end;
      $reserved     = 0; 
      $block_count  = 1;
      $block_sizes  = $$featR{$uid}->{length}.",";
      $block_starts = "0,";
      @bedline = join ("\t", ($chrom,$chrom_start,$chrom_end,
			      $name,$score,$strand,$thick_start,
			      $thick_end,$reserved,$block_count,
			      $block_sizes, $block_starts));
      print BEDOUT "@bedline\n";
    }
   close (BEDOUT);

   # sort bed file
   my $cmd = "$bedtools sort -i $bedname_u > $bedname";
   my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
     run( command => $cmd, verbose => 0 );
   if( !$success ) {
     print STDERR "ERROR [$this_function] Call to $bedtools  unsuccessful\n";
     print STDERR "ERROR: this is what the command printed:\n";
     print join "", @$full_buf;
     croak $!;
   }
   unlink($bedname_u);

 } # end foreach
 if (defined $log){close(LOG)};
}

# feature_summary($fsR,$dest)
# Print summary of %featstat hash
#
# ARG1: reference to %featstat hash
# ARG2: path for output file
sub feature_summary {
  my ($fsR, $dest) = @_;
  my ($fn,$fh);
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] $dest does not exist\n"
    unless (-d $dest);

  $fn = dir($dest,$$fsR{accession}.".summary.txt");
  open $fh, ">", $fn or croak $!;

  print $fh "Accession\t $$fsR{accession}\n";
  print $fh "Origin   \t $$fsR{origin}\n";
  foreach my $ft (sort keys %$fsR){
    next if ($ft =~ /total/ || $ft =~ /accession/ || $ft =~ /origin/);
    print $fh "$ft\t$$fsR{$ft}\n";
  }
  print $fh "Total\t$$fsR{total}\n";
  close $fh;
}


1;

__END__

=head1 NAME

Bio::ViennaNGS::AnnoC - Perl extension for converting sequence
annotation formats

=head1 SYNOPSIS

  use Bio::ViennaNGS::AnnoC;

  parse_gff($gff3_file);

  feature_summary($fstat,$dest);

  features2bed($featR,$fstatR,$feature,$dest,$bn,$log)

  get_fasta_ids($fasta_file);

=head1 DESCRIPTION

=over 4

=item parse_gff($gff3_file)

C<parse_gff()> parses GFF3 annotation files of non-spliced
genomes. The GFF3 specification is available at
L<http://www.sequenceontology.org/resources/gff3.html> This routine
expects the path to a GFF3 file as argument C<$gff3_file> and returns
two hash references: C<$feat> is a reference to a HOH containing the
raw annotation information for each feature found in the GFF3
file. C<$fstat> references a hash containing summary statistics of the
features found in the GFF3 file.

This routine has been tested with NCBI bacteria GFF3 annotation.

=item feature_summary($fstat,$dest)

This routine generates a summary file for all features parsed by
parse_gff. It expects two arguments: C<$fstat> is a refence to the
summary hash generated by C<parse_gff()> and C<$dest> is the output
path for a summary.txt file.

=item features2bed($featR,$fstatR,$feature,$dest,$bn,$log)

Thsi routine converts genomic features from a L<Bio::ViennaNGS::AnnoC>
C<%features> hash into BED12 format. C<$featR> and C<$fstatR> are
references to the C<%features> and C<%featstat> hashes,
respectively. C<$feature> can be either a string corresponding to a
genbank key present in C<%features> or C<undef>. If it is defined,
only features of the speficied key will be written to one single BED12
file. If C<$feature> is undef, BED12 files will be generated for each
genbank key present in C<%features>. C<$dest> is the output directory
and C<$bn> the basename for all output files. C<$log> can either be
the full path to a logfile or C<undef>.

=item get_fasta_ids($fasta_file)

C<get_fasta_ids()> returns an array containing all headers/IDs of a
Fasta file as L<Bio::DB:Fasta> objects. The Fasta file may contain
multiple entries.

=back

=head1 DEPENDENCIES

=over 4

=item L<Bio::Tools::GFF>

=item L<Bio::DB::Fasta>

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
