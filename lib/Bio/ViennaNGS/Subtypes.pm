# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 19:16:16 michl>

package Bio::ViennaNGS::Subtypes;

use version; our $VERSION = qv('0.17');
use Moose::Util::TypeConstraints;
use Bio::DB::Fasta;
use Params::Coerce ();

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
