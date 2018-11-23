# -*-CPerl-*-
# Last changed Time-stamp: <2018-11-23 17:43:57 mtw>

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
  message { "[SUBTYPE ISSUE]: expecting '+','-', or '.', however ($_) is neither."};

subtype 'Bio::ViennaNGS::ZeroOrOne',
  as 'Int',
  where { $_ == 0 or $_ == 1 },
  message { "[SUBTYPE ISSUE]: expecting '0' or '1', however ($_) is neither." };
