#include <Sequence/polySiteVector8.hpp>
#include <Sequence/PolyTable.hpp>
#include <iostream>
#include <algorithm>
namespace Sequence
{
  polySiteVector8::polySiteVector8( const std::initializer_list<polymorphicSite> & p) : base()
  {
    std::cerr << "polySiteVector8 ilist\n";
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first, Seq8(p.second,dna_poly_alphabet)) );
		   });
  }
  
  polySiteVector8::polySiteVector8( const PolyTable & p) : base()
  {
    std::cerr << "polySiteVector8 const PolyTable\n";
    this->reserve(p.numsites());
    std::for_each( p.sbegin(), p.send(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first, Seq8(p.second,dna_poly_alphabet)) );
		     //poly8::dna2vtype(p.second)) );
		   });
  }

  polySiteVector8::polySiteVector8( const std::vector<polymorphicSite> & p) : base()
  {
    std::cerr << "polySiteVector8 const vector psite\n";
    this->reserve(p.size());
    std::for_each( p.cbegin(), p.cend(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first, Seq8(p.second,dna_poly_alphabet)) );//poly8::dna2vtype(p.second)) );
		   });
  }
  
  polySiteVector8::polySiteVector8( std::vector<polymorphicSite> & p) : base()
  {
    std::cerr << "polySiteVector8 vector psite\n";
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( polymorphicSite & p ) {
		     this->push_back( polymorphicSite8(p.first,Seq8(p.second,std::ref(dna_poly_alphabet))) );//poly8::dna2vtype(p.second)) );
						//poly8::dna2vtype(std::ref(p.second))) );
		   });
    p.clear();
  }
}
