# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-05 15:25:05 michl>

package Bio::ViennaNGS::Subtypes;

use version; our $VERSION = qv('0.17_03');
use Moose::Util::TypeConstraints;
use Bio::DB::Fasta;
use Params::Coerce ();


subtype 'Bio::ViennaNGS::MyFasta' => as class_type('Bio::DB::Fasta');

coerce 'Bio::ViennaNGS::MyFasta'
  => from 'Str'
  => via { Bio::DB::Fasta->new($_) }
  => from 'Object'
  => via {$_ -> isa('Bio::DB::Fasta') ? $_ : Params::Coerce::coerce ('Bio::DB::Fasta', $_); }

