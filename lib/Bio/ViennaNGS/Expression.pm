# -*-CPerl-*-
# Last changed Time-stamp: <2015-02-11 16:52:47 mtw>

package Bio::ViennaNGS::Expression;

use version; our $VERSION = qv('0.12');
use Moose;
use Carp;
use Data::Dumper;
use Path::Class;
use Bio::ViennaNGS::Bed;
use Bio::ViennaNGS::Util qw(sortbed);

use namespace::autoclean;


has 'readcountfile' => (
		    is => 'rw',
		    predicate => 'has_readcountfile',
		  );

has 'data' => (
	       is => 'rw',
	       isa => 'ArrayRef',
	       default => sub { [] },
		   );

has 'conds' => (
		is => 'rw',
		isa => 'Int',
		predicate => 'has_conds',
	       );

has 'nr_features' => (
		      is => 'rw',
		      isa => 'Int',
		      predicate => 'has_features',
		     );


sub parse_readcounts_bed12 {
  my ($self,$file) = @_;
  my @mcData = ();
  my ($i,$n) = 0x2;
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] readcount / multicov file $self->readcountfile not available\n"
    unless (-e $file);
  $self->readcountfile($file);
  open (RC_IN, "< $file") or croak $!;

  while (<RC_IN>){
    $n++;
    chomp;
    # 0:chr|1:start|2:end|3:name|4:score|5:strand
    # 6:thickStart|7:thickEnd|8:itemRgb|9:blockCount|
    # 10:blockSizes|11:blockStarts
    @mcData = split(/\t/);
    my $conditions = (scalar @mcData)-12;  # multicov extends BED12
    $self->conds($conditions);

    # NOTE: Better keep BED12 entries in a hash, generating UUIDs as
    # keys instead of storing the same BED12 entry n times (ie for
    # each sample) in $self->data

    my $bedobj =  Bio::ViennaNGS::Bed->new(chromosome   => $mcData[0],
					   start        => $mcData[1],
					   end          => $mcData[2],
					   name         => $mcData[3],
					   score        => $mcData[4],
					   strand       => $mcData[5],
					   thickStart   => $mcData[6],
					   thickEnd     => $mcData[7],
					   itemRgb      => $mcData[8],
					   blockCount   => $mcData[9],
					   blockSizes   => $mcData[10],
					   blockStarts  => $mcData[11],
					  );
    my $len = $bedobj->length;
    my $id = sprintf("%s:%d-%d_%s_%s", $mcData[0],$mcData[1],
		     $mcData[2],$mcData[3],$mcData[5]);
    #print "\$id: $id\n";
    for ($i=0;$i<$self->conds;$i++){
      ${$self->data}[$i]{$id} = {
				 bed_entry => $bedobj,
				 length    => $len,
				 count     => $mcData[eval(12+$i)],
				};
      # print Dumper(${$self->data}[$i]);
    }
  }
  $self->nr_features($n);
  close(RC_IN);
}

sub write_expression_bed12 {
  my ($self,$item,$dest,$base_name) = @_;
  my ($bedname,$bedname_u,$outfile,$feat,$i);
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function]: $dest does not exist\n"
    unless (-d $dest);

  $bedname = $base_name.".".$item.".multicov.bed12";
  $bedname_u = $base_name.".".$item.".multicov.u.bed12";

  $outfile = file($dest,$bedname_u);

  croak "ERROR [$this_function]: $self->conds not available\n"
    unless ($self->has_conds);

  croak "ERROR [$this_function]: $self->nr_features not available\n"
    unless ($self->has_features);

  #print "=====> write_multicov: writing multicov file $bedfile with ".
  #  eval($self->nr_features)."  lines and ".eval($self->conds)." conditions\n";

  # check whether each element in @{$self->data} has the same amount of entries
  for($i=0;$i<$self->conds;$i++){
    my $fc = scalar keys %{$self->data->[$i]}; # of keys in %{$featCount}
    #print "condition $i => $fc keys\n";
    croak "ERROR [$this_function]: unequal element count in @{$self->data}"
      unless($self->nr_features == $fc)
  }

  open (MULTICOV_OUT, "> $bedname_u") or croak $!;

  # use BED12 data stored with condition 0 here, assuming its the same for all conditions
  foreach $feat (keys %{${$self->data}[0]} ){
    my @mcLine = ();
    # retrieve BED12 line first
    my $bedo = ${${$self->data}[0]}{$feat}{bed_entry};
    my $bedline = $bedo->as_bed_line(12);
    push @mcLine, $bedline;

    # process multicov values for all samples
    for($i=0;$i<$self->conds;$i++){
      # print "------------>  ";  print "processing $i th condition ";  print "<-----------\n";
      croak "Could not find item $feat in mcSample $i\n"
	unless ( defined ${${$self->data}[$i]}{$feat} );
      push @mcLine, ${${$self->data}[$i]}{$feat}{$item};

    }
   print MULTICOV_OUT join("\t",@mcLine)."\n";
  }
  close(MULTICOV_OUT);

  sortbed($bedname_u,"./",$bedname,1,undef);  # sort bed file
}


sub computeTPM {
  my ($self,$sample,$rl) = @_;
  my ($TPM,$T,$totalTPM) = (0)x3;
  my ($i,$meanTPM);

  # iterate through $self->data[$i] twice:
  # 1. for computing T (denominator in TPM formula)
  foreach $i (keys %{${$self->data}[$sample]}){
    my $count  = ${${$self->data}[$sample]}{$i}{count};
    my $length =  ${${$self->data}[$sample]}{$i}{length};
    #print "count: $count\nlength: $length\n";
    $T += $count * $rl / $length;
  }

  # 2. for computng actual TPM values
  foreach $i (keys %{${$self->data}[$sample]}){
    my $count  = ${${$self->data}[$sample]}{$i}{count};
    my $length =  ${${$self->data}[$sample]}{$i}{length};
    $TPM = 1000000 * $count * $rl/($length * $T);
    ${${$self->data}[$sample]}{$i}{TPM} = $TPM;
    $totalTPM += $TPM;
  }

  $meanTPM = $totalTPM/$self->nr_features;

 # print Dumper(${$self->data}[$sample]);
 # print "totalTPM=$totalTPM | meanTPM=$meanTPM\n";


  return $meanTPM;
}


__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::Expression - An object oriented interface for
read-count based gene expression

=head1 SYNOPSIS

  use Bio::ViennaNGS::Expression;

  my $expression = Bio::ViennaNGS::Expression->new();

  # parse read counts from an extended BED12 file
  $expression->parse_readcounts_bed12("$bed12");

  # compute normalized expression of ith sample in Transcript per Million (TPM)
  $expression->computeTPM($i, $readlength);

  # write extended BED12 file with TPM for each condition past
  # the 12th column
  $expression->write_expression_bed12("TPM", $dest, $basename);

=head1 DESCRIPTION

This module provides a L<Moose> interface for computation of gene /
transcript expression from read counts.

=head1 METHODS

=over

=item parse_readcounts_bed12

Title : parse_readcounts_bed12

Usage : C<$obj-E<gt>parse_readcounts_bed12($file);>

Function : Parses a bedtools multicov (multiBamCov) file, i.e. an
           extended BED12 file, into an Array of Hash of Hashes data
           structure (C<@{$self-E<gt>data}>).

Args : C<$file> is the input file, i.e. and extended BED12 file where
       each column past the 12th lists read counts for this bedline's
       feature(s) for a specific sample/condition.

Returns :

Notes : This method evaluates the number of samples/conditions present
        in the input, i.e. the number of columns extending the
        canonical BED12 columns in the input multicov file and
        populates C<$self-E<gt>conds>. Also populates
        C<$self-E<gt>nr_features> with the number of genes/features
        present in the input (evidently, this should be the same for
        each sample/condition in the input).

=item computeTPM

Title : computeTPM

Usage : C<$obj-E<gt>computeTPM($sample, $readlength);>

Function : Computes expression of each gene/feature present in
           C<$self-E<gt>data> in Transcript per Million (TPM) [Wagner
           et.al. Theory Biosci. (2012)].  is a reference to a Hash of
           Hashes data straucture where keys are feature names and
           values hold a hash that must at least contain length and
           raw read counts. Practically, C<$featCount_sample> is
           represented by _one_ element of C<@featCount>, which is
           populated from a multicov file by C<parse_multicov()>.

Args : C<$sample> is the sample index of C<@{$self-E<gt>data}>. This is
        especially handy if one is only interested in computing
        normalized expression values for a specific sample, rather
        than all samples in multicov BED12 file. C<$readlength> is the
        read length of the RNA-seq sequencing experiment.

Returns : Returns the mean TPM of the processed sample, which is
          invariant among samples. (TPM models relative molar
          concentration and thus fulfills the invariant average
          criterion.)

=item write_expression_bed12

Title : write_expression_bed12

Usage : C<$obj-E<gt>write_expression_bed12($measure,
$dest,$basename);>

Function : Writes normalized expression data to a bedtools multicov
           (multiBamCov)-type BED12 file.

Args : C<$measure> specifies the type in which normalized expression
       data from C<@{$self-E<gt>data}> is dumped. Currently supports
       'TPM', however 'RPKM' support will be added in a future
       release. These values must have been computed and inserted into
       C<@{self-E<gt>data}> beforehand by
       e.g. C<$self-E<gt>computeTPM()>. C<$dest> and C<$base_name>
       give path and base name of the output file, respectively.

Returns : None. The output is position-sorted extended BED12 file.

=back

=head1 DEPENDENCIES

=over

=item L<Moose>

=item L<Carp>

=item L<Path::Class>

=item L<namespace::autoclean>

=back

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::Bed>

=item L<Bio::ViennaNGS::Util>

=back

=head1 AUTHOR

Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
