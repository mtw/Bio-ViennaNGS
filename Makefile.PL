use 5.010000;
use lib 'inc';
use ExtUtils::MakeMaker;
use strict;
use warnings;
use File::ShareDir::Install;

install_share dist => 'share';

my @scripts = qw (assembly_hub_constructor.pl bam_split.pl bam_to_bigwig.pl bam_uniq.pl bed2bedGraph.pl bed2nt2aa.pl bed62bed12.pl gff2bed.pl fasta_subgrep.pl fasta_multigrep.pl kmer_analysis.pl newUCSCdb.pl normalize_multicov.pl rnaseq_peakfinder.pl sj_visualizer.pl splice_site_summary.pl track_hub_constructor.pl trim_fastq.pl Tutorial_pipeline03.pl);
WriteMakefile(
	      'NAME'            => 'Bio::ViennaNGS',
	      'AUTHOR'          => q{Michael T. Wolfinger},
	      'VERSION_FROM'    => 'lib/Bio/ViennaNGS.pm',
	      'ABSTRACT_FROM'   => 'lib/Bio/ViennaNGS.pm',
	      'LICENSE'         => 'perl',
	      'PREREQ_PM'       => {
				    "Getopt::Long"	           => '2.0',
				    "constant"                     => '0',
				    "Carp" 		           => '0',
				    "Cwd"                          => '0',
				    "Bio::Perl"                    => '0',
				    "Bio::Root::Version"           => '1.00690001',
				    "Bio::SeqUtils"                => '0',
				    "Bio::DB::Sam"                 => '1.37',
				    "Bio::DB::Fasta"               => '0',
				    "Bio::Tools::GFF"              => '0',
				    "Data::Dumper"                 => '0',
				    "Env"                          => '0',
				    "Exporter"                     => '0',
				    "File::Basename"     	   => '0',
				    "File::Find"                   => '0',
				    "File::Path"                   => '0',
				    "File::Share"                  => '0',
				    "File::Temp" 	           => '0',
				    "File::Slurp"                  => '0',
				    "FindBin"                      => '0',
				    "IPC::Cmd" 		           => '0',
				    "lib"                          => '0',
				    "List::Util"                   => '0',
				    "Moose"      		   => '0',
				    "MooseX::Clone"                => '0',
				    "Moose::Util::TypeConstraints" => '0',
				    "namespace::autoclean"         => '0',
				    "Path::Class" 	           => '0',
				    "Math::Round"                  => '0',
				    "Params::Coerce"               => '0',
				    "Pod::Usage" 	           => '0',
				    "PerlIO::gzip"                 => '0',
				    "POSIX" 		           => '0',
				    "strict"                       => '0',
				    "Template"                     => '0',
				    "Test::Files"                  => '0',
				    "Test::File::Contents"         => '0',
				    "Test::Deep"                   => '0',
				    "Tie::Hash::Indexed"           => '0',
				    "version"                      => '0',
				    "warnings"                     => '0',
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
