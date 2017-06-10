# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-10 19:13:43 michl>

package Bio::ViennaNGS::Subtypes;

use Bio::ViennaNGS;
use Moose::Util::TypeConstraints;
use Bio::DB::Fasta;
use Params::Coerce ();
use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");

subtype 'Bio::ViennaNGS::MyFasta' => as class_type('Bio::DB::Fasta');

coerce 'Bio::ViennaNGS::MyFasta'
  => from 'Str'
  => via { Bio::DB::Fasta->new($_) }
  => from 'Object'
  => via {$_ -> isa('Bio::DB::Fasta') ? $_ : Params::Coerce::coerce ('Bio::DB::Fasta', $_); };

subtype 'Bio::ViennaNGS::PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- nor ."};
