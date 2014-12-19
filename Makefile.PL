use 5.012004;
use lib 'inc';
use ExtUtils::MakeMaker;
use strict;
use warnings;
use File::ShareDir::Install;

install_share dist => 'share';

my @scripts = qw (assembly_hub_constructor.pl bam_quality_stat.pl bam_split.pl bam_to_bigwig.pl bam_uniq.pl bed2bedGraph.pl extend_bed.pl gff2bed.pl kmer_analysis.pl MEME_XML_motif_extractor.pl motiffinda.pl newUCSCdb.pl normalize_multicov.pl sj_visualizer.pl splice_site_summary.pl trim_fastq.pl Tutorial_pipeline01.pl);
WriteMakefile(
	      'NAME'            => 'Bio::ViennaNGS',
	      'AUTHOR'          => q{Michael T. Wolfinger},
	      'VERSION_FROM'    => 'lib/Bio/ViennaNGS.pm',
	      'ABSTRACT_FROM'   => 'lib/Bio/ViennaNGS.pm',
	      'LICENSE'         => 'perl',
	      'PREREQ_PM'       => {
				    "Getopt::Long"	        => '2.0',
				    "Bio::Root::Version" 	=> '1.00690001',
				    "Bio::DB::Sam"              => '1.39',
				    "Bio::DB::Fasta"            => '0',
				    "Bio::Tools::GFF"           => '0',
				    "Data::Dumper"              => '0',
				    "File::Basename"     	=> '0',
				    "File::Temp" 	        => '0',
				    "IPC::Cmd" 		        => '0',
				    "Path::Class" 	        => '0',
				    "Math::Round"               => '0',
				    "Pod::Usage" 	        => '0',
				    "PerlIO::gzip"              => '0',
				    "POSIX" 		        => '0',
				    "Carp" 		        => '0',
				    "Tie::Hash::Indexed"        => '0',
				    "Moose"      		=> '0',
				    "MooseX::Clone"             => '0',
				    "Moose::Util::TypeConstraints" => 0,
				    "namespace::autoclean"  => '0',
				    "File::Share"           => '0',
				   },
	      'MIN_PERL_VERSION'          => '5.010',
	      'EXE_FILES'                 => [ map { "scripts/$_" } @scripts ],
	      'CONFIGURE_REQUIRES'        => {
					      'ExtUtils::MakeMaker' => 0,
					     },
	      'BUILD_REQUIRES' => {
				   'Test::More' => 0,
				  },
	      dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
	      clean => { FILES => 'Bio-ViennaNGS-*' },
	     );

package MY;
use File::ShareDir::Install qw(postamble);
