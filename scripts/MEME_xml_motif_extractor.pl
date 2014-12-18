#!/usr/bin/env perl
#Last changed Time-stamp: <2014-12-16 12:26:58 fall>
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
use Storable;
use Cwd;
use All::Misc;
use File::Path qw(make_path remove_tree);
use XML::Simple;
use Statistics::R;
###############
###Variables
###############

my $VERBOSE = 0;
my ( $dir, $odir, $file, $outname, $discard, $meme, $totalsites, $acid, $rlibpath );
my ( @matches );
my (%set, %motifs, %regex, %sequences, %seq );
###############
###Command Line Options
###############

pod2usage(-verbose => 0)
	unless GetOptions(
	"dir|d=s"     => \$dir,
	"odir|o=s"    => \$odir,
	"file|f=s"    => \$file,
	"outname|t=s" => \$outname,
	"discard|x=s" => \$discard,
	"acid|a=s"    => \$acid,
	"rlib|r=s"    => \$rlibpath,
	"help|h"      => sub{pod2usage(-verbose => 1)},
	"man|m"       => sub{pod2usage(-verbose => 2)},
	"verbose"     => sub{ $VERBOSE++ }
	);

$dir	 =  cwd() unless ($dir);
$odir	 =  "$dir"."/XMLParseOut" unless $odir;
$dir	 =~ s/ //g;
$odir	 =~ s/ //g;
$outname =  "" unless (defined $outname);
$outname =~ s/\.[bed|fa]*//g if ($outname);
$acid    = 'DNA' unless (defined $acid && $acid eq 'RNA');
($dir) or die "No working directory chosen!\n";

die 'No R-library path defined! Find out using the comman \'.libPaths()\' from an R shell!\n' unless ($rlibpath);
die 'No xml file defined! Please provide a valid meme.xml file with the -f option!\n' unless ($file);
if (defined $discard){
    @matches=split(/,/,$discard);
}
else{
    @matches = qw();
}

##############
###  Main  ###
##############

if (!-d $odir){
make_path($odir) or die "Error creating directory: $odir";
}

chdir($dir) or die "Directory $dir not found!\n";

$meme = XMLin($file);
%seq = %{$meme->{training_set}->{sequence}};

foreach my $sequence ( keys %seq ){
    my $seq_id = $seq{"$sequence"}{id};
    $set{$seq_id} = $sequence;
    my $name = $sequence;
    $sequences{$seq_id}= $name;
    $totalsites++;
}

my %startpos;
my %motif = %{$meme->{motifs}};
foreach my $mo ( keys %{$motif{motif}} ){
    (my $regtmp = $motif{motif}{$mo}{regular_expression}) =~ s/\n//g;
    $regex{$mo} = $regtmp;
    $regex{$mo} =~ s/T/U/g if ($acid eq 'RNA');
    my $regx = $regex{$mo};
    my @sites = @{$motif{motif}{$mo}{contributing_sites}{contributing_site}};
    my $evalue = $motif{motif}{$mo}{e_value};
    my $width = $motif{motif}{$mo}{width};

    for (0..$#sites){
	my %entry = %{$sites[$_]};
	my $seqid = $entry{sequence_id};
	my $seq = $sequences{$seqid};
	my $re = qr/$regx/;
	while($seq =~ /$re/g){
	    push @{$startpos{$seqid}} , pos($seq);
	}
	my $start = $entry{position};
	my $end   = $start + $width -1;
	push @{$motifs{$regx}{seq}}, $seqid;
	$motifs{$regx}{sites}++;
	$motifs{$regx}{nr}=$mo;
	$motifs{$regx}{ev}=$evalue;
	push @{$motifs{$regx}{start}}, $start;
	push @{$motifs{$regx}{end}}, $end;
    }
}

chdir($odir) or die "Directory $odir not found!\n";

open (OUT, ">", "Motif_Site_Distribution$outname");
print OUT "Motif\tRegex\tSites\tEvalue\n";
print OUT "Total\tTotal\t$totalsites\t\n";

my @mof = (sort {$motifs{$b}{sites} <=> $motifs{$a}{sites}} keys %motifs );
foreach my $motf (@mof){
    my $sites = $motifs{$motf}{sites};
    my $ex = $motf;
    my $nr = $motifs{$motf}{nr};
    my $ev = $motifs{$motf}{ev};
    print OUT "Motif$nr\t$ex\t$sites\t$ev\n";
}
close (OUT);

print STDERR "Running R to plot site distribution!\n";
my $R = Statistics::R->new();
$R->startR;
$R->set('rlib', $rlibpath);
$R->run(q`.libPaths( c( .libPaths(), rlib) )`);
$R->run(q`x<-.libPaths()`);
$R->run(q`library(ggplot2);`);
$R->run(q`library(RColorBrewer);`);
$R->set('file', "Motif_Site_Distribution".$outname);
$R->run(q`dat <- read.table(file, head=T, sep="\t",colClasses=c('character','character','numeric','character'))`);
$R->run(q`data <- head(dat,n=11)`);
$R->run(q`motif<-data$Motif`);
$R->run(q`regex<-data$Regex`);
$R->run(q`site<-data$Sites`);
$R->run(q`eval<-data$Evalue`);
$R->run(q`colourSig = length(unique(regex))`);
$R->run(q`getPalette = colorRampPalette(brewer.pal(30, "Set1"))`);
$R->run(q`p <- ggplot(data, aes(x=reorder(regex, site), y=site, fill=regex),guide=FALSE) + geom_bar(stat = "identity", position = "stack", guide=FALSE) + scale_fill_manual(values = getPalette(colourSig),guide=FALSE) + coord_flip() + geom_text(aes(label=eval), position=position_dodge(width=2), size=4, hjust=0) + scale_y_continuous(expand=c(0.15,0))`);
$R->run(q`p <- p + theme(aspect.ratio=1)`);
$R->run(q`p <- p + theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold"))`);
$R->run(q`p <- p + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))`);
$R->run(q`p <- p + theme(axis.title.y = element_text(angle=90))`);
$R->run(q`p <- p + ylab("Number of sites")`);
$R->run(q`p <- p + xlab("Motif regular expression")`);
$R->run(q`p <- p + ggtitle(file)+theme(plot.title = element_text(size=12, face="bold"))`);
$R->run(q`p`);
$R->run(q`out<-paste(file,".eps",sep="")`);
$R->run(q`ggsave(filename=out, path="./", width=12, height=8)`);
$R->run(q`out<-paste(file,".svg",sep="")`);
$R->run(q`ggsave(filename=out, path="./", width=12, height=8)`);
$R->run(q`out<-paste(file,".jpg",sep="")`);
$R->run(q`ggsave(filename=out, path="./", width=12, height=8)`);
$R->stopR;

if(-e "Rplots.pdf"){ ### Clean outdir
    `rm Rplots.pdf`;
}

chdir($dir) or die "Directory $dir not found!\n";

##################################END################################

###############
###POD
###############

__END__

=head1 NAME

MEME_xml_motif_extractor.pl - Generates a simple statistic with ggplots about the sequence coverage of motifs from MEME standard xml output

=head1 SYNOPSIS
MEME_xml_motif_extractor.pl [-f I<FILE>] [-d I<STRING>] [-o I<STRING>] [-t I<STRING>] [-x I<String>] [-a I<STRING>] [-r I<STRING>]

=head1 OPTIONS

=over 

=item B<-f>

MEME xml output file

=item B<-d>

Path to working dir containing the xml file

=item B<-o>

Path to Output directory

=item B<-t>

Name of output files

=item B<-x>

If annotation is found in the sequence name fields, items that should not be analyzed (e.g. annnotion like miscRNA, ncRNA, predicted_gene) can be defined here

=item B<-a>

Defined the type of nucleic acid in MEME input, defaults to DNA, can be ignored if protein sequences were used

=item B<-r>

Path to the R library

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 DESCRIPTION

This program digests MEME xml output and returns a list of found motifs with the number of sequences containing those motifs

=head1 AUTHOR

Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

=cut


##################################END################################


