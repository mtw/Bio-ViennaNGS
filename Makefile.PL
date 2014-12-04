use 5.012004;
use ExtUtils::MakeMaker;
my @scripts = qw (bam_split.pl bam_to_bigwig.pl bam_uniq.pl gff2bed.pl motiffinda.pl newUCSCdb.pl normalize_multicov.pl sj_visualizer.pl splice_site_summary.pl);
WriteMakefile(
    NAME              => 'Bio::ViennaNGS',
    VERSION_FROM      => 'lib/Bio/ViennaNGS.pm',
    PREREQ_PM         => {
			  "Getopt::Long" => 2.0,
			  "Bio::Perl" => 1.6,
			  "Bio::DB::Sam"  => 1.39,
			  "Data::Dumper" => 0,
			  "File::Basename" => 0,
			  "File::Temp" => 0,
			  "IPC::Cmd" => 0,
			  "Path::Class" => 0,
			  "Pod::Usage" => 0,
			  "POSIX" => 0,
			  "Carp" => 0,
			  "Moose" => 0,
			 },
  EXE_FILES           => [ map { "scripts/$_" } @scripts ],
    ($] >= 5.005 ?
      (ABSTRACT_FROM  => 'lib/Bio/ViennaNGS.pm',
       AUTHOR         => 'Michael T. Wolfinger <michael@wolfinger.eu>') : ()),
);
