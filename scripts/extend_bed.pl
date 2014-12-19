#!/usr/bin/env perl
# Last changed Time-stamp: <2014-12-19 23:45:09 mtw>
# AUTHOR: Joerg Fallmann <joerg.fallmann@univie.ac.at>

###############
###Use stuff
###############
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use PerlIO::gzip;
use Cwd;
use File::Path qw(make_path remove_tree);
use File::Basename qw(fileparse);
use Math::Round;
use Bio::ViennaNGS::Util qw(extend_chain parse_bed6);
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::FeatureChain;

###############
###Variables
###############

my $VERBOSE = 0;
my ( $g, $b, $o, $l, $r, $e, $u, $d );

###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
    "genome|g=s"  => \$g,
    "bedfile|b=s" => \$b,
    "outfile|o=s" => \$o,
    "left|l=s"    => \$l,
    "right|r=s"   => \$r,
    "extend|e=s"  => \$e,
    "up|u=s"      => \$u,
    "down|d=s"    => \$d,
    "help|h"      => sub{pod2usage(-verbose => 1)},
    "man|m"       => sub{pod2usage(-verbose => 2)},      
    "verbose"     => sub{ $VERBOSE++ }
    );

###############
### MAIN
###############

unless (-f $g && -f $b ) {
	warn "Could not open bed or genome file for processing, please provide them via the -g and -f options!\n";
	pod2usage(-verbose => 0);
}

unless ( $e || $l || $r || $u || $d){ 
    warn "No number of flanking nucleotides chosen, output would be input! Please provide them via the -l, -r, -e, -u or -d shortoption\n";
    pod2usage(-verbose => 0);
}

if ($e){
    $l=nearest(1,$e/2);
    $r=nearest(1,$e/2);
}

$l=0 unless $l;
$r=0 unless $r;
$d=0 unless $d;
$u=0 unless $u;

open (my $Genome, "<:gzip(autopop)",$g) or die "$!";

my $filextension;
if ($d){
    $filextension .= "_".$d."_fromEnd";
}
if ($u){
    $filextension .= "_".$u."_fromStart";
}
if ($l){
    $filextension .= "_".$l."_upstream";
}
if ($r){
    $filextension .= "_".$r."_downstream";
}
if ($e){
    $filextension = "_".$e."_equal";
}

my($filename, $dirs, $suffix) = fileparse($b,'.bed');
print STDERR "fil:$filename\tdir:$dirs\tsuf:$suffix\n";
$o = $filename.$filextension.$suffix;

open (my $Out, ">",$o) or die "$!";

my %sizes;

while (<$Genome>){
    chomp $_;
    my ($chr,$size)=split (/\t/,$_);
    $sizes{$chr}=$size;
}

my @featurelist = @{parse_bed6($b)};
### Now we create a Bio::ViennaNGS::FeatureChain from the featurelist above
my $chain = Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);

### This chain is now processed, and all features are extended and stored into a separate chain
my $extended_chain = extend_chain(\%sizes,$chain,$l,$r,$u,$d);

print STDERR "Creating Output";
my $out = $extended_chain->print_chain();
print $Out $out;

###############
###POD
###############

__END__

=head1 NAME

extend_bed.pl - Extends bed entries strand specific one- or two-sided.

=head1 SYNOPSIS

extend_bed.pl [-g I<FILE>] [-b I<FILE>] [-o I<FILE>] [-e I<Interger>]
[-l I<Interger>] [-r I<Interger>] [-u I<Interger>] [-d I<Interger>]
[options]

=head1 DESCRIPTION

This program extends Bed files to a total length of at least -e, -r or
-l nucleotides, or at least to begin or end of the corresponding
chromosome. Furthermore 3' or 5' end extension is possible, where only
up- or downstream regions can be retrieved, or combinations of the
latter with extension into the original coordinate span.

=head1 OPTIONS

=over 

=item B<-g>

File containing chromosome sizes

=item B<-b>

BED file for extension

=item B<-o>

Output file name

=item B<-e>

Extension to total length from both sides

=item B<-l>

Extension to total length from 5' end

=item B<-r>

Extension to total length from 3' end

=item B<-u>

Upstream only extension to total length from 5' end. Can be combined
with -d to get upstream extension, plus extension into original
coordinate span

=item B<-d>

Downstream only extension to total length from 3' end.  Can be
combined with -l to get downstream extension, plus extension into
original coordinate span

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################

