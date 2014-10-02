use 5.012004;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Bio::ViennaNGS',
    VERSION_FROM      => 'lib/Bio/ViennaNGS.pm', # finds $VERSION
    PREREQ_PM         => {
			  "Bio::Perl" => 0,
			  "Bio::DB::Sam"  => 1.39,
			  "Data::Dumper" => 0,
			  "File::Basename" => 0,
			  "File::Temp" => 0,
			  "IPC::Cmd" => 0,
			  "Path::Class" => 0,
			  "Carp" => 0,
			 }, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Bio/ViennaNGS.pm', # retrieve abstract from module
       AUTHOR         => 'Michael T. Wolfinger <michael@wolfinger.eu>') : ()),
);
