#!/usr/bin/env perl
# Last changed Time-stamp: <2014-12-10 12:18:26 mtw>

###############
###Use stuff
###############
use strict;
use warnings;
use PerlIO::gzip;
use Cwd;
use File::Path qw (make_path);
use Pod::Usage;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Bio::ViennaNGS qw(unique_array);

###############
###Variables
###############

my $VERBOSE=0;
my ( $file, $type, $peak, $chromsize, $track, $conv, $anno );

###############
###Command Line Options
###############

pod2usage(-verbose => 0)
    unless GetOptions(
	"file|f=s"	=> \$file,
	"type|t=s"	=> \$type,
	"chrom|c=s"	=> \$chromsize,
	"anno|a=s"	=> \$anno,
	"help|h"	=> sub{pod2usage(-verbose => 1)},
	"man|m"		=> sub{pod2usage(-verbose => 2)},
	"verbose"	=> sub{ $VERBOSE++ }
    );


my $dir = cwd();

my $file = file("Bedgraph", $file);

make_path($file) unless (-d $file);

die "You must provide filename and cell-type, and file with chrom-sizes\n"
  unless ($file && $type && $chromsize);

###############
###MAIN
###############

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";
print STDERR "Adding annotation information\n" if ($anno && $anno eq "extended");
my %sizes=();

open (SIZE, "<:gzip(autopop)", $chromsize);
while (<SIZE>){
    chomp (my $raw = $_);
    my @line = split("\t",$raw);
    $sizes{"$line[0]"}=$line[1];
}
close (SIZE);

open (IN, "<:gzip(autopop)" , $file);
chdir("Bedgraph/$file") or die "Could not move to Bedgraph $!\n";

my %covplus  = ();
my %covminus = ();
my %anno     = ();
my %annop    = ();
my %annom    = ();

while (<IN>){
    chomp (my $raw = $_);
    push my @line , split (/\t/,$raw);
    push @line, "\." if ( !$line[5] ); 
    my $chrom  = $line[0];#)=~ s/chr//g;
    my $start  = $line[1];
    my $end    = $line[2];
    my $strand = $line[5];
    my $annotation = join("\t",$line[3],$line[4]);
    for (5..$#line){
	$annotation .= "\t".$line[$_];
    }
    for (my $i=$start;$i<=$end;$i++){
	my $a=$i+1;
	next if ($a > ($sizes{"$chrom"}));
	if ($strand eq "+" or $strand eq "1"){
	    $covplus{"$chrom\t$i\t$a"}+=1;
	    push @{$annop{"$chrom\t$i\t$a"}}, $annotation if ($anno && $anno eq "extended");
	}
	elsif($strand eq "-" or $strand eq "-1"){
	    $covminus{"$chrom\t$i\t$a"}+=1;
	    push @{$annom{"$chrom\t$i\t$a"}}, $annotation if ($anno && $anno eq "extended");
	}
    }
}
close (IN);

foreach my $key (keys %covplus){
    my $cov=$covplus{$key};
    my @an = unique_array(\@{$annop{$key}});
    my $annotation = join("\t",@an);
    my @tmp=split(/\t/,$key);
    my $chrom=$tmp[0];
    open (OUT1, ">>:gzip", "$chrom\.$type\.fw.gz");
    if ($anno && $anno eq "extended"){
	print OUT1 "$key\t$cov\t$annotation\n";
    }
    else{
	print OUT1 "$key\t$cov\n";
    }
    close (OUT1);
}

foreach my $key (keys %covminus){
    my $cov=$covminus{$key}; # at least one occurence
    my @an = unique_array(\@{$annom{$key}});
    my $annotation = join("\t",@an);
    my @tmp=split(/\t/,$key);
    my $chrom=$tmp[0];
    open (OUT2, ">>:gzip", "$chrom\.$type\.re.gz");
    if ($anno && $anno eq "extended"){
	print OUT2 "$key\t-$cov\t$annotation\n";
    }
    else{
	print OUT2 "$key\t-$cov\n";
    }
    close (OUT2);
}

`for i in \`ls chr\*\.fw.gz\`;do zcat \$i|sort -k1,1 -k2,2 -n|gzip > \$i\.bedg.gz;done`;
`for i in \`ls chr\*\.re.gz\`;do zcat \$i|sort -k1,1 -k2,2 -n|gzip > \$i\.bedg.gz;done`;

chdir ($dir) or die "Could not move to $dir $!\n";


###############
###POD
###############

__END__

=head1 NAME

bed2bedgraph.pl - Convert BED or extended BED files to
bedGraph format

=head1 SYNOPSIS

bed2bedgraph.pl [-f I<FILE>] [-c I<FILE>] [-t I<STRING>] [-a
I<STRING>] [options]

=head1 DESCRIPTION

This program converts BED files to strand specific bedGraph files,
allowing additional annotation and automatic generation of bedGraph
files which can easily be converted to big-type files for UCSC
visualization.

=head1 OPTIONS

=over 

=item B<-f>

BED file for conversion.

=item B<-c>

File containing chromosome sizes

=item B<-t>

Type of bed file (e.g. sample name, replicate name, condition, ...)

=item B<-a>

Whether file is a standard bed or extended bed, 'extended' for extended bed

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back



=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################

