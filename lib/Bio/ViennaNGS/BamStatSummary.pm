# -*-CPerl-*-
# Last changed Time-stamp: <2015-01-05 15:49:29 mtw>

package Bio::ViennaNGS::BamStatSummary;

use 5.12.0;
use version; our $VERSION = qv('0.12_08');
use Moose;
use Carp;
use POSIX qw(floor);
use Statistics::R;
use Data::Dumper;
use Path::Class;
use namespace::autoclean;
use Bio::ViennaNGS::BamStat;
use Tie::Hash::Indexed;
use File::Basename;

has 'data' => (
	       is => 'ro',
	       isa => 'ArrayRef [Bio::ViennaNGS::BamStat]',
	       default => sub { [] },
	      );

has 'countStat' => (
		    is => 'rw',
		    isa => 'HashRef',
		    predicate => 'has_countStat',
		    default => sub { {} },
#		    auto_deref => '1',
		   );

has 'outpath' => (
		  is => 'rw',
		  isa => 'Str',
		  required => '1',
		 );

has 'rlib' => (
	       is => 'rw',
	       isa => 'Str',
	       required => '1',
	      );

has 'files' => (
		is => 'rw',
		isa => 'ArrayRef',
		required => 1,
		predicate => 'has_files',
	       );

has 'control_match' => ( # provides stats how many mapped bases match the reference genome
			is => 'rw',
			isa => 'Bool',
			default => '1',
			predicate => 'has_control_match',
		       );

has 'control_clip' => ( # provides stats how many bases are soft or hard clipped
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       predicate => 'has_control_clip',
		       );

has 'control_split' => ( # provides stats how many/often  mapped reads are split
			is => 'rw',
			isa => 'Bool',
			default => '1',
			predicate => 'has_control_split',
		       );

has 'control_qual' => ( # provides stats on quality of the match
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       predicate => 'has_control_qual',
		      );

has 'control_edit' => ( # provides stats on the edit distance between read and mapped reference
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       predicate => 'has_control_edit',
		      );

has 'control_flag' => ( # analyses the sam bit flag for qual/strands/pair_vs_single reads
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       predicate => 'has_control_flag',
		      );

has 'control_score' => ( # provides stats on per-base quality scores
			is => 'rw',
			isa => 'Bool',
			default => '1',
			predicate => 'has_control_score',
		       );

has 'control_uniq' => ( #  gives number and stats of multiplicity of readaligments
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       predicate => 'has_control_uniq',
		      );

has 'is_segemehl' => ( # toggles to consider segemehl specific bam feature
		      is => 'rw',
		      isa => 'Bool',
		      default => '0',
		      predicate => 'has_is_segemehl',
		     );

sub populate_data {
  my ($self) = @_;
  foreach my $bamfile (@{$self->files}){
    carp ">> processing $bamfile\n";
    my $bo = Bio::ViennaNGS::BamStat->new(bam => $bamfile);
    $bo->stat_singleBam();
    push (@{$self->data}, $bo);
  }
}


sub populate_countStat {
  my ($self) = @_;

  tie my %hdr, 'Tie::Hash::Indexed';
  %hdr = (
	  "sample"           => "Sample",
	  "total_alignments" => "# Total alignments",
	  "mapped_reads"     => "# Mapped reads",
	  "umapped_reads"    => "# Unique mapped reads",
	  "mmapped_reads"    => "# Multi mapped reads",
	  "aligned_pairs"    => "# Aligned in pairs",
	  "aligned_mm"       => "# Aligned mate missing",
	  "aligned_se"       => "# Aligned single end",
	  "aligned_fwd"      => "# Aligned forward strand",
	  "aligned_rev"      => "# Aligned reverse strand");
  ${$self->countStat}{'header'} = \%hdr;

  foreach my $sample (@{$self->data}){
    my ($basename,$dir,$ext) = fileparse($$sample{'bam'},qr/\.[^.]*/);
     ${$self->countStat}{$basename}{'total_alignments'} = floor(0.5 + $$sample{'data_out'}->{'aln_count'}->{'total'} );
     ${$self->countStat}{$basename}{'mapped_reads'} = floor(0.5 + $$sample{'data_out'}->{'uniq'}->{'mapped_reads'} );
     ${$self->countStat}{$basename}{'umapped_reads'} = floor(0.5 + $$sample{'data_out'}->{'uniq'}->{'uniq_mapped_reads'} );
     ${$self->countStat}{$basename}{'mmapped_reads'} = 
       floor(0.5 + $$sample{'data_out'}->{'uniq'}->{'mapped_reads'} - $$sample{'data_out'}->{'uniq'}->{'uniq_mapped_reads'} );
     ${$self->countStat}{$basename}{'aligned_pairs'} = floor(0.5 + ($$sample{'data_out'}->{'aln_count'}->{'mapped_pair'})/2 );
     ${$self->countStat}{$basename}{'aligned_mm'} = floor(0.5 + $$sample{'data_out'}->{'aln_count'}->{'unmapped_pair'} );
     ${$self->countStat}{$basename}{'aligned_se'} = floor(0.5 + $$sample{'data_out'}->{'aln_count'}->{'mapped_single'} );
     ${$self->countStat}{$basename}{'aligned_fwd'} = floor(0.5 + $$sample{'data_out'}->{'strand'}->{'forward'} );
     ${$self->countStat}{$basename}{'aligned_rev'} = floor(0.5 + $$sample{'data_out'}->{'strand'}->{'reverse'} );
  }
}

sub dump_countStat {
  my ($self,$how) = @_;
  my $mn = "mapping_stats.csv";
  my $fn = file($self->outpath,$mn);
  open(OUT, "> $fn") or croak "cannot open OUT $!";

  print OUT join ("\t", values %{$self->countStat->{'header'} })."\n";
  foreach my $sample (keys %{$self->countStat} ){
    next if ($sample eq 'header');
    my @line = ();
    foreach my $key ( keys %{$self->countStat->{'header'}} ) {
      if ($key eq 'sample'){
	push (@line, $sample);
	next;
      }
      push @line, $self->countStat->{$sample}->{$key};
    }
    print OUT join ("\t", @line)."\n";
  }

#  print Dumper($self->countStat);
  close (OUT);
}

sub make_BarPlot{
  my ($self) = @_;
  my @Rstat_data_count = ();

  ## collect data for read.table string
  push @Rstat_data_count, 'Samples', grep {!/header/} keys %{$self->countStat}; # first line with sample names
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                              # end 1st line

  push @Rstat_data_count, 'aligned_se';
  foreach my $sample ( grep {!/header/}  keys %{$self->countStat} ) {    # 2nd line with single end aln counts
    push @Rstat_data_count, $self->countStat->{$sample}{'aligned_se'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                       # end 2nd line

  push @Rstat_data_count, 'aligned_pairs';
  foreach my $sample ( grep {!/header/}  keys %{$self->countStat} ) {    # 3rd line with aligned pairs count
    push @Rstat_data_count, $self->countStat->{$sample}{'aligned_pairs'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                       # end 3rd line

  push @Rstat_data_count, 'aligned_mm';
  foreach my $sample ( grep {!/header/}  keys %{$self->countStat} ) {    # 4th line with incomplete aligned pairs
    push @Rstat_data_count, $self->countStat->{$sample}{'aligned_mm'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                       # end 4th line

  ## produce bar plot
  my $mn = "mapping_stats.pdf";
  my $fn = file($self->outpath,$mn);
  my $datastring = join(" ", @Rstat_data_count);
  $self->plot_barplot($fn, "Mapped reads", $datastring); # produce plot with read.table string input
}

sub plot_barplot { #plot barplot read.table text string
  my ($self, $filename,$ylab,$data_string) = @_;
  my ($bn,$odir,$ext) = fileparse($filename, qr /\..*/);
  #my $rlibpath        = '/usr/bin/R';
  my $rlibpath        = $self->rlib;

  $filename .= '.pdf' unless ($ext eq 'pdf');

  my $R = Statistics::R->new();
  $R->startR;
  $R->set('rlib', $rlibpath);
  $R->set('log_dir', $odir);
  $R->run("pdf('${filename}')") ;
  $R->run("dat<-read.table(text = \"$data_string\", header = TRUE, row.names=1)") ;
  $R->run("dat_m<-as.matrix(dat)") ;
  $R->run("colors<-terrain.colors(nrow(dat_m), alpha = 1)") ;
# $R->run("colors<-c('darkolivegreen4','darkolivegreen2','coral1')") ;
# $R->run("colors<-hcl(seq(0, 360, length =nrow(dat_m)))") ;
  $R->run("barplot(dat_m, xlim=c(0,ncol(dat_m)+2), col=colors, legend.text = TRUE, args.legend = list(x = ncol(dat_m) + 2, y=max(colSums(dat_m)), bty = 'n' ), ylab='$ylab', xlab='Samples')") ;
  $R->run("dev.off()") ;
  $R->stopR;
}

1;
