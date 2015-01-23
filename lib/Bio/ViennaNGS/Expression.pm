# -*-CPerl-*-
# Last changed Time-stamp: <2015-01-22 15:54:54 mtw>

package Bio::ViennaNGS::Expression;

use version; our $VERSION = qv('0.12_11');
use Moose;
use Carp;
use Data::Dumper;

use namespace::autoclean;


has 'readcountfile' => (
		    is => 'ro',
		    isa => 'Str',
		    
		  );

has 'data' => (
	       is => 'rw',
	       isa => 'ArrayRef',
	       predicate => 'has_data',
		   );

has 'conds' => (
		is => 'rw',
		isa => 'Int',
		lazy => 1,
	       );


sub parse_readcounts_bed12 {
  my ($self) = @_;
  my @mcData = ();
  my $i;
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function] readcount / multicov file $self->readcountfile not available\n"
    unless (-e $file);
  open (RC_IN, "< $file") or croak $!;

  while (<RC_IN>){
    chomp;
    @mcData = split(/\t/); # 0:chr|1:start|2:end|3:name|4:score|5:strand
    $self->conds = (scalar @mcData)-6; # multicov extends BED6
    #print "$_\n";
    for ($i=0;$i<$self->conds;$i++){
      ${$self->data}[$i]{$mcData[3]} = {
				    chr    => $mcData[0],
				    start  => $mcData[1],
				    end    => $mcData[2],
				    name   => $mcData[3],
				    score  => $mcData[4],
				    strand => $mcData[5],
				    len    => eval($mcData[2]-$mcData[1]),
				    count  => $mcData[eval(6+$i)],
				   }
    }
    #print Dumper(@mcData);
  }
  close(RC_IN);
  return $mcSamples;
}

sub write_multicov {
  my ($item,$dest,$base_name) = @_;
  my ($outfile,$mcSamples,$nrFeatures,$feat,$i);
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function]: $dest does not exist\n"
    unless (-d $dest);
  $outfile = file($dest,$base_name.".".$item.".multicov.csv");
  open (MULTICOV_OUT, "> $outfile") or croak $!;

  $mcSamples = scalar @featCount; # of samples in %{$featCount}
  $nrFeatures = scalar keys %{$featCount[1]}; # of keys in %{$featCount}[1]
  #print "=====> write_multicov: writing multicov file $outfile with $nrFeatures lines and $mcSamples conditions\n";

  # check whether each column in %$featCount has the same number of entries
  for($i=0;$i<$mcSamples;$i++){
    my $fc = scalar keys %{$featCount[$i]}; # of keys in %{$featCount}
    #print "condition $i => $fc keys\n";
    unless($nrFeatures == $fc){
      croak "ERROR [$this_function]: unequal element count in \%\$featCount\nExpected $nrFeatures have $fc in condition $i\n";
    }
  }

  foreach $feat (keys  %{$featCount[1]}){
    my @mcLine = ();
    # process standard BED6 fields first
    push @mcLine, (${$featCount[1]}{$feat}->{chr},
		   ${$featCount[1]}{$feat}->{start},
		   ${$featCount[1]}{$feat}->{end},
		   ${$featCount[1]}{$feat}->{name},
		   ${$featCount[1]}{$feat}->{score},
		   ${$featCount[1]}{$feat}->{strand});
    # process multicov values for all samples

    for($i=0;$i<$mcSamples;$i++){
     # print "------------>  ";  print "processing $i th condition ";  print "<-----------\n";
      unless (defined ${$featCount[$i]}{$feat}){
	croak "Could not find item $feat in mcSample $i\n";
      }
      push @mcLine, ${$featCount[$i]}{$feat}->{$item};

    }
    #print Dumper(\@mcLine);
    print MULTICOV_OUT join("\t",@mcLine)."\n";
  }
  close(MULTICOV_OUT);
}


sub computeTPM {
  my ($featCount_sample,$rl) = @_;
  my ($TPM,$T,$totalTPM) = (0)x3;
  my ($i,$features,$meanTPM);

  $features = keys %$featCount_sample; # of of features in hash
  print Dumper(\%$featCount_sample);print ">>$rl<<\n";

  # iterate through $featCount_sample twice:
  # 1. for computing T (denominator in TPM formula)
  foreach $i (keys %$featCount_sample){
    $T += ($$featCount_sample{$i}{count} * $rl)/($$featCount_sample{$i}{len});
  }
  # 2. for computng actual TPM values
  foreach $i (keys %$featCount_sample){
    $TPM = 1000000 * $$featCount_sample{$i}{count} * $rl/($$featCount_sample{$i}{len} * $T);
    $$featCount_sample{$i}{TPM} = $TPM;
    $totalTPM += $TPM;
  }
  $meanTPM = $totalTPM/$features;
  # print "totalTPM=$totalTPM | meanTPM=$meanTPM\n";
  return $meanTPM;
}


__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::Expression - OO interface for read-count based gene
expression

=head1 SYNOPSIS

  use Bio::ViennaNGS::Expression;

  my $file = 

=head1 DESCRIPTION

This module provides a L<Moose> interface for computation of gene /
transcript expression based on read counts.

=head1 METHODS

=head1 DEPENDENCIES

=over 3

=item L<Moose>

=item L<Carp>

=item L<namespace::autoclean>

=back

=head1 SEE ALSO

=over 1

=item L<Bio::ViennaNGS>

=back

=head1 AUTHOR

Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
