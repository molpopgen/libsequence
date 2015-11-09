#include <Sequence/polySiteVector8.hpp>
#include <Sequence/PolyTable.hpp>
#include <algorithm>

namespace Sequence
{
  polySiteVector8::polySiteVector8( const std::initializer_list<polymorphicSite> & p) : base()
  {
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first, Seq8(p.second,dna_poly_alphabet)) );
		   });
  }
  
  polySiteVector8::polySiteVector8( const PolyTable & p) : base()
  {
    this->reserve(p.numsites());
    std::for_each( p.sbegin(), p.send(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first, Seq8(p.second,dna_poly_alphabet)) );
		   });
  }

  polySiteVector8::polySiteVector8( const std::vector<polymorphicSite> & p) : base()
  {
    this->reserve(p.size());
    std::for_each( p.cbegin(), p.cend(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first, Seq8(p.second,dna_poly_alphabet)) );
		   });
  }
}
