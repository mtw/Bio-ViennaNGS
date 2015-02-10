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
my ($track_hub_return_value,$lf);
my $logname = "log.txt";
my $genome_identifier = 'hg19';
my $folder_in = '-';
my $dest = '.';
my $base_URL = '-';
my $chrom_size_file = '-';
my $big_wig_urls = '-';
###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
    "infolder|i=s" => \$folder_in,
    "out|o=s"      => \$dest,
    "baseurl|b=s"  => \$base_URL,
    "chromsize|c=s"  => \$chrom_size_file,
    "bigwigs|bw=s" => \$big_wig_urls,
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

  perl Tutorial_Pipeline03.pl [--infolder I<PATH>] [--out I<PATH>] [--baseurl -I<URL>] [--bigwigs -I<URL,URL#URL>]

  This tutorial is based on the track_hub_constructor.pl script and the output from Tutorial02_pipeline.pl. While the option descriptions here are specific for the results from tutorial02 the trackhub_hub_constructor.pl can be applied in the same manner to other datasets.
  Example call Tutorial_Pipeline03.pl --infolder /scratch/egg/projects/viennangs/tutorial02_output --out /scratch/egg/projects/viennangs/tutorial03 --baseurl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/ --bigwigs http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw#http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw

=head1 DESCRIPTION

This script demonstrates UCSC genome browser trackhub construction with <Bio::ViennaNGS>.

The results of Tutalorial02_pipeline.pl are vizualized, showing genomic region, annotation
and expression level, which can be interpreted in conjunction to each other. 

=head2 PREREQUITES

For running this tutorial on your machine you will need a full
installation of the L<Bio::ViennaNGS> distribution, and result files
obtained from tutorial02.pl, which can also be downloaded from 
L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS>):

=over

=item F<hg19_highlyexpressed.bed >

=item F<hg19_highlyexpressed.ext50_upstream.bed >

=item F<hg19_highlyexpressed.pos.bed>

=item F<hg19_highlyexpressed.neg.bed>

=item F<hg19_highlyexpressed.pos.bw>

=item F<hg19_highlyexpressed.neg.bed>

=item F<hg19.chrom.sizes>

=back

For the UCSC genome browser to be able to
use our trackhub it needs to be accessible via URL. This is the base
URL you need to provide as commandline argument. 

=head2 DISCLAIMER

The resulting trackhub can only be read by the UCSC genome browser if it is accessible via URL.
If you have no webspace for testing available you can use example output available from
our webserver L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS>.


=head1 PIPELINE

=head3 Create UCSC Genome Browser Trackhub

=cut

print STDERR "Constructing UCSC genome browser trackhub ...\n";

$track_hub_return_value = make_track_hub($genome_identifier,$folder_in,$dest,$base_URL,$chrom_size_file,$big_wig_urls,$lf);


print STDERR "DONE\n";

=head1 COMMAND LINE OPTIONS

=cut

__END__


=over 4

=head1 OPTIONS

=over


=item B<--infolder -i>

Directory which contains all track files in BED/bigBed format. The
resulting Track Hub will contain these files in their respective
bigFile version.

=item B<--out -o>

Destination folder for the output Track Hub.

=item  B<--baseurl -b>

BaseURL used within the Track Hub. This URL will be included verbatim
in the resulting Track Hub. It is crucial that this URl is valid, else
the resulting Track Hub will be broken.

=item  B<--bigwigs -bw>

URLs pointing to big wig files to be included in the trackhub. Multiple URLs are
separated by the star character #. It is possible to create a multiwig container by
providing 2 URLs instead of one separated by comma character ,. E.g.
http://foo.com/bar.bw,http://foo.com/bar2.bw#http://foo.com/bar3.bw yields a multi
big wig container displaying bar as positive reads in green and bar2 as negative
3 red colored reads in the same track and additionally bar3 in an own track
colored blue.

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
