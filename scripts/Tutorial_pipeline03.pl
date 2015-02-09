#!/usr/bin/env perl
# Last changed Time-stamp: <2015-02-09 17:02:46 egg>
# AUTHOR: Florian Eggenhofer <florian.eggenhofer@univie.ac.at>

###############
###Use stuff
###############
use strict;
use warnings;
use Cwd;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use File::Path qw(make_path remove_tree);
use Bio::ViennaNGS::UCSC qw( make_assembly_hub make_track_hub );
use IPC::Cmd qw(can_run run);
###############
###Variables
###############

my $VERBOSE = 0;
my $base_URL = '-';
###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
    "baseurl|b=s"  => \$base_URL,
    "help|h"    => sub{pod2usage(-verbose => 1)},
    "man|m"     => sub{pod2usage(-verbose => 2)},
    "verbose"   => sub{ $VERBOSE++ }
    );

###############
### MAIN
###############

=head1 NAME

Tutorial_Pipeline03.pl - Construct a UCSC genome browser trackhub

=head1 SYNOPSIS

  perl Tutorial_Pipeline03.pl [--infolder I<PATH>] [--out
I<PATH>] [--baseurl -I<URL>] [--bigwigs -I<URL,URL#URL>]
  Where base_url URL will be included verbatim in the resulting Track Hub. 
  It is crucial that this URl is valid, else the resulting Track Hub will be broken.
  If no base_url is provided, trackhub creation will be skipped.

=head1 DESCRIPTION

This script demonstrates UCSC genome browser trackhub construction with <Bio::ViennaNGS>.



=head2 PREREQUITES

For running this tutorial on your machine you will need a full
installation of the L<Bio::ViennaNGS> distribution, and result files
obtained from tutorial02.pl, which can also be downloaded from 
L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS>):

=over

=item F<hg19_chromchecked.fa>

=item F<hg19_highlyexpressed.bed>

=item F<hg19.chrom.sizes>

=item F<C1R1.bam>

=item F<hg19.chrom.sizes>

=back

=head2 DISCLAIMER

The resulting trackhub can only be read by the UCSC genome browser if it is accessible via URL.
If you have no webspace for testing available you can use example output available from
our webserver L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS>.


=head1 PIPELINE

Let's first initialize some variables and generate a chromosome_sizes
hash.
  my ($track_hub_return_value,$lf);
  my $logname = "log.txt";
  my $genome_identifier = '-';
  my $folder_in = '-';
  my $dest = '.';
  my $base_URL = '-';
  my $chrom_size_file = '-';
  my $big_wig_urls = '-';

=cut

  my ($track_hub_return_value,$lf);
  my $logname = "log.txt";
  my $genome_identifier = '-';
  my $folder_in = '-';
  my $dest = '.';
  my $base_URL = '-';
  my $chrom_size_file = '-';
  my $big_wig_urls = '-';

=head3 Create UCSC Genome Browser Trackhub

Construction of a UCSC Genome Browser Trackhub trackhub vizualizes the
results. Genomic region, annotation and expression level can be interpreted
in conjunction to each other. For the UCSC genome browser to be able to
use our trackhub it needs to be accessible via URL. This is the base
URL you need to provide as commandline argument. 

=cut

print STDERR "Constructing UCSC genome browser trackhub ...\n";

$track_hub_return_value = make_track_hub($genome_identifier,$folder_in,$dest,$base_URL,$chrom_size_file,$big_wig_urls,$lf);


print STDERR "DONE\n";

=head1 COMMAND LINE OPTIONS

=cut

__END__


=over 4

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHORS

=over

=item Florian Eggenhofer E<lt>florian.eggenhofer@univie.ac.atE<gt>

=back

=cut


##################################END################################
