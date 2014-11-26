#!/bin/perl

#Script ExtendBed.pl;
#Last changed Time-stamp: <2014-11-26 15:09:17 fall> by Joerg Fallmann <joerg.fallmann@univie.ac.at>

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
use Math::Round;
use Bio::ViennaNGSutil qw(extend_chain parse_bed6);
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::FeatureIO;
###############
###Variables
###############

my $VERBOSE = 0;
my ( $g, $b, $o, $l, $r, $e );

###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
	"genome|g=s"     => \$g,
	"bedfile|b=s"    => \$b,
	"outfile|o=s"    => \$o,
	"left|l=s"       => \$l,
	"right|r=s"      => \$r,
	"extend|e=s"     => \$e,
	"help|h"         => sub{pod2usage(-verbose => 1)},
	"man|m"          => sub{pod2usage(-verbose => 2)},      
	"verbose"        => sub{ $VERBOSE++ }
    );

###############
###MAIN
###############

unless (-f $g && -f $b ) {
	warn "Could not open bed or genome file for processing, please provide them via the -g and -f options!\n";
	pod2usage(-verbose => 0);
}

unless ( $e || $l || $r ){ 
    warn "No number of flanking nucleotides chosen, output would be input!Please provde them via the -l and -r or -e option\n";
    pod2usage(-verbose => 0);
}

if ($e){
    $l=nearest(1,$e/2);
    $r=nearest(1,$e/2);
}
$l=0 unless $l;
$r=0 unless $r;

print STDERR "Using $g to extend $b by $l and $r and save it to $o\n";

open (my $Genome, "<:gzip(autopop)",$g) or die "$!";
#open (my $Bed, "<:gzip(autopop)",$b) or die "$!";
open (my $Out, ">",$o) or die "$!";

my %sizes;

while (<$Genome>){
    chomp $_;
    my ($chr,$size)=split (/\t/,$_);
    $sizes{$chr}=$size;
}
### This depends on MTW, either FeatureIO reads files and distributes to objects, or not
### So far I do it directly in this script
### Possible workarounds would be 
#### my $bedobject = FeatureIO->new($b,);
#### $bedobject->featurechain_from_Bed($b); ## Reads bed file into featurechain object

#my @featurelist; ## This will become a FeatureChain
#while(<$Bed>){
##    ### This should be done by FeatureIO if want to;
#    chomp (my $raw = $_);
#    push my @line , split (/\t/,$raw);
#    push @line, "\." if ( !$line[5] ); 
#
#    (my $chromosome  = $line[0])=~ s/chr//g;
#    my $start	     = $line[1]+1;
#    my $end	     = $line[2];
#    my $name	     = $line[3];
#    my $score	     = $line[4];
#    my $strand	     = $line[5];
#    my $extension = '';
#
#    if ($line[6]){
#	for (6..$#line){
#	    $extension .= $line[$_]."\t";
#	}
#	$extension = substr($extension,0,-1);
#    }
#    my $feat = Bio::ViennaNGS::Feature->new(chromosome=>$chromosome,start=>$start,end=>$end,name=>$name,score=>$score,strand=>$strand,extension=>$extension);
#    push @featurelist, $feat;
#}
#
my @featurelist = @{parse_bed6($b)};
### Now I create a Bio::ViennaNGS::FeatureChain from the featurelist above
my $chain = Bio::ViennaNGS::FeatureChain->new('type'=>'original','chain'=>\@featurelist);
#print STDERR "Dumping Chain of length ".scalar($#{$chain->chain})."\n";
#print Dumper($chain);

### This chain is now processed, and all features are extended and stored into a separate chain
print STDERR "Extending with $l and $r\n";
my $extended_chain = extend_chain(\%sizes,$chain,$l,$r);
#print STDERR "Dumping ExChain  of length ".scalar($#{$chain->chain})."\n";
#print Dumper($extended_chain);

print STDERR "Creating Output";
my $out = $extended_chain->print_chain();
print $Out $out;

###############
###POD
###############

__END__

=head1 NAME

ExtendBed.pl - Extends bed entries strand specific one- or two-sided.

=head1 SYNOPSIS
ExtendBed.pl [-g I<FILE>] [-b I<FILE>] [-o I<FILE>] [-e I<Interger>] [-l I<Interger>] [-r I<Interger>]
[options]

=head1 OPTIONS

=over 

=item B<-g>

File containing chromosome sizes

=item B<-b>

Bed file for extension

=item B<-o>

Output file name

=item B<-e>

Extension to total length from both sides

=item B<-l>

Extension to total length from 5' end

=item B<-r>

Extension to total length from 3' end

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 DESCRIPTION

This program extends Bed files to a total length of at least -e, -r or -l nucleotides, or at least to begin or end of the corresponding chromosome

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################

