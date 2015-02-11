#!/usr/bin/env perl
# Last changed Time-stamp: <2015-02-11 14:27:50 fall>
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
use File::Basename;
use Path::Class;
use Bio::ViennaNGS::UCSC qw( make_track_hub );
###############
###Variables
###############

my $VERBOSE = 0;
my ($track_hub_return_value,$lf);
my $logname = "log.txt";
my $genome_identifier = 'hg19';
my $dest = '.';
my $base_URL = '-';
my $big_bed_urls = '';
my $big_wig_urls = '';
###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
    "out|o=s"      => \$dest,
    "baseurl|u=s"  => \$base_URL,
    "bigbeds|b=s" => \$big_bed_urls,
    "bigwigs|w=s" => \$big_wig_urls,
    "help|h"    => sub{pod2usage(-verbose => 1)},
    "man|m"     => sub{pod2usage(-verbose => 2)},
    "verbose"   => sub{ $VERBOSE++ }
    );

###############
### MAIN
###############
print "$base_URL\n";
unless ($base_URL =~ /^http/) {
  warn "Base URL must be given in full format, eg http://foo.bar.com/Hubs/";
  pod2usage(-verbose => 0);
}

unless ($dest =~ /\/$/){$dest.= "/";}
unless (-d $dest){
  mkdir $dest or die $!;
}
$lf = file($dest,$logname);

=head1 NAME

Tutorial_Pipeline03.pl - Construct a UCSC genome browser trackhub

=head1 SYNOPSIS

  perl Tutorial_Pipeline03.pl [--out I<PATH>] [--baseurl -I<URL>] [--bigbeds -I<URL#URL>] [--bigwigs -I<URL,URL#URL>]

=head1 DESCRIPTION

This script demonstrates UCSC genome browser trackhub construction with <Bio::ViennaNGS>.

The results of Tutalorial02_pipeline.pl are vizualized in UCSC,
showing genomic region, annotation and expression level, which can be
interpreted in conjunction to each other.

The result of this tutorial can be viewed by navigating your browser to 
L<here|http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hubUrl=http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_trackHub/trackHub/hub.txt&position=chr15>

This tutorial is based on the track_hub_constructor.pl script and the output from Tutorial02_pipeline.pl.
While the option descriptions here are specific for the results from Tutorial02 the trackhub_hub_constructor.pl
can be applied in the same manner to other datasets. The example call uses the bigwig and bedfiles available from our server.
Example call: Tutorial_pipeline03.pl -o /home/user/public_html/hg19_trackHub -u http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_trackHub -b http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bb#http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bb -w http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw,http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw

=head2 PREREQUITES

For running this tutorial on your machine you will need a full
installation of the L<Bio::ViennaNGS> distribution, and result files
obtained from tutorial02.pl, which can alternatively be downloaded
from L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS>):

=over

=item F<hg19_highlyexpressed.pos.bb>

=item F<hg19_highlyexpressed.neg.bb>

=item F<hg19_highlyexpressed.pos.bw>

=item F<hg19_highlyexpressed.neg.bw>

=back

For the UCSC genome browser to be able to
visualize our trackhub it needs to be accessible via URL. This is the base
URL you need to provide as commandline argument. 

=head2 DISCLAIMER

The resulting trackhub can only be read by the UCSC genome browser if it is accessible via URL.
If you have no webspace for testing available you can use example output available from
our webserver L<here|http://nibiru.tbi.univie.ac.at/ViennaNGS>.


=head1 PIPELINE

=head3 Create UCSC Genome Browser Trackhub

=cut

print "Constructing UCSC genome browser trackhub ...\n";

$track_hub_return_value = make_track_hub($genome_identifier,$dest,$base_URL,$big_bed_urls,$big_wig_urls,$lf);

print "Loading the trackhub into the UCSC genome browser...\n\n";

print "1. Point your browser to: http://genome.ucsc.edu/index.html\n\n";

print "2. Select Genome Browser from the left menu\n\n";

print "3. You are now redirected to nearest mirror of the genome browser.\n    Select \"My Data\" from the top menu and then \"Track Hubs\" in the popup list.\n\n";

print "4. The \"Track Data Hubs\" page is displayed. In the register \"My Hubs\" it is possible to add the newly created track hub. Paste the URL to hub.txt into the URL field and click \"Add Hub\" (e.g. http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_trackHub/trackHub/hub.txt).\n\n";

print "5. The track hub now loads into the hg19 public hub. \n\n";

print "6. Enter e.g. chr15 in the position field and hit go\n\n";

print "7. You should now see 2 annotation tracks and a bigwig multi-track plotted in red and green\n\n";

print "DONE\n";

=head1 COMMAND LINE OPTIONS

=cut

__END__


=over 4

=head1 OPTIONS

=over

=item B<--out -o>

Destination folder for the output Track Hub.

=item  B<--baseurl -u>

BaseURL used within the Track Hub. This URL will be included verbatim
in the resulting Track Hub. It is crucial that this URl is valid, else
the resulting Track Hub will be broken.

=item  B<--bigbeds -b>

URLs pointing to big bed files to be included in the trackhub. Multiple URLs are
separated by the character #. 

=item  B<--bigwigs -w>

URLs pointing to big wig files to be included in the trackhub. Multiple URLs are
separated by the character #. It is possible to create a multiwig container by
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
