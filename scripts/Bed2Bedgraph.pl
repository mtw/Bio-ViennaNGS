#!/usr/bin/perl
# Last changed Time-stamp: <2014-11-06 16:01:49 fall>

##############
###Library for Testing, remove before delivery
##############
use lib '/scratch/fall/Work/ViennaNGS/viennangs/lib';

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
use Bio::ViennaNGSutil qw(unique_array);

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
	"peak|p=s"	=> \$peak,
	"chrom|c=s"	=> \$chromsize,
	"track|a=s"	=> \$track,
	"converted|v=s" => \$conv,
	"anno|o=s"	=> \$anno,
	"help|h"	=> sub{pod2usage(-verbose => 1)},
	"man|m"		=> sub{pod2usage(-verbose => 2)},      
	"verbose"	=> sub{ $VERBOSE++ }
    );


my $dir = cwd();

if (!-d "Bedgraph/$file"){
    make_path("Bedgraph/$file");
}

die "You must provide filename and cell-type, and file with chrom-sizes\n" unless ($file && $type && $chromsize);

###############
###MAIN
###############

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";
print STDERR "Analyzing Peak bedfile\n" if ($peak || $file =~ /bed_spk/);
print STDERR "Generating track\n" if ($track && $track eq "track");
print STDERR "Analyzing already genomic positions\n" if ($conv && $conv eq "on");
print STDERR "Adding annotation information\n" if ($anno && $anno eq "anno");
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
	    push @{$annop{"$chrom\t$i\t$a"}}, $annotation if ($anno && $anno eq "anno");
	}
	elsif($strand eq "-" or $strand eq "-1"){
	    $covminus{"$chrom\t$i\t$a"}+=1;
	    push @{$annom{"$chrom\t$i\t$a"}}, $annotation if ($anno && $anno eq "anno");
	}
    }
}
close (IN);

foreach my $key (keys %covplus){
    if ($track && $track eq "track"){
	my $cov=$covplus{$key};
	my @an = unique_array(\@{$annop{$key}});
	my $annotation = join("\t",@an);
	my @tmp=split(/\t/,$key);
	my $chrom=$tmp[0];
	open (OUT1, ">>:gzip", "$chrom\.$type\.fw.track.gz");
	if ($anno && $anno eq "anno"){
	    print OUT1 "$key\t$cov\t$annotation\n";
	}
	else{
	    print OUT1 "$key\t$cov\n";
	}
	close (OUT1);
    }
    else{
	my $cov=$covplus{$key};
	my @an = unique_array(\@{$annop{$key}});
	my $annotation = join("\t",@an);
	my @tmp=split(/\t/,$key);
	my $chrom=$tmp[0];
	open (OUT1, ">>:gzip", "$chrom\.$type\.fw.gz");
	if ($anno && $anno eq "anno"){
	    print OUT1 "$key\t$cov\t$annotation\n";
	}
	else{
	    print OUT1 "$key\t$cov\n";
	}
	close (OUT1);
    }
}

foreach my $key (keys %covminus){
    if ($track && $track eq "track"){
	my $cov=$covminus{$key};
	my @an = unique_array(\@{$annom{$key}});
	my $annotation = join("\t",@an);
	my @tmp=split(/\t/,$key);
	my $chrom=$tmp[0];
	open (OUT1, ">>:gzip", "$chrom\.$type\.re.track.gz");
	if ($anno && $anno eq "anno"){
	    print OUT1 "$key\t-$cov\t$annotation\n";
	}
	else{
	    print OUT1 "$key\t-$cov\n";
	}
	close (OUT1);
    }
    else{
	my $cov=$covminus{$key}; # at least one occurence
	my @an = unique_array(\@{$annom{$key}});
	my $annotation = join("\t",@an);
	my @tmp=split(/\t/,$key);
	my $chrom=$tmp[0];
	open (OUT2, ">>:gzip", "$chrom\.$type\.re.gz");
	if ($anno && $anno eq "anno"){
	    print OUT2 "$key\t-$cov\t$annotation\n";
	}
	else{
	    print OUT2 "$key\t-$cov\n";
	}
	close (OUT2);
    }
}

if ($track && $track eq "track"){
    `for i in \`ls chr\*\.fw.track.gz\`;do zcat \$i|sort -k1,1 -k2,2 -n|gzip > \$i\.bedg.track.gz;done`;
    `for i in \`ls chr\*\.re.track.gz\`;do zcat \$i|sort -k1,1 -k2,2 -n|gzip > \$i\.bedg.track.gz;done`;
}
else{
    `for i in \`ls chr\*\.fw.gz\`;do zcat \$i|sort -k1,1 -k2,2 -n|gzip > \$i\.bedg.gz;done`;
    `for i in \`ls chr\*\.re.gz\`;do zcat \$i|sort -k1,1 -k2,2 -n|gzip > \$i\.bedg.gz;done`;
}
chdir ($dir) or die "Could not move to $dir $!\n";
