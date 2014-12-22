#include <Sequence/Ptable8.hpp>
#include <Sequence/PolyTable.hpp>
#include <iostream>
namespace Sequence
{
  Ptable8::Ptable8( const std::initializer_list<polymorphicSite> & p) : base()
  {
    std::cerr << "Ptable8 ilist\n";
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, Seq8(p.second,dna_poly_alphabet)) );
		   });
  }
  
  Ptable8::Ptable8( const PolyTable & p) : base()
  {
    std::cerr << "Ptable8 const PolyTable\n";
    this->reserve(p.numsites());
    std::for_each( p.sbegin(), p.send(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, Seq8(p.second,dna_poly_alphabet)) );
		     //poly8::dna2vtype(p.second)) );
		   });
  }

  Ptable8::Ptable8( const PolyTable * p) : base()
  {
    std::cerr << "Ptable8 const PolyTable *\n";
    *this = Ptable8(*p);
  }

  Ptable8::Ptable8( const std::vector<polymorphicSite> & p) : base()
  {
    std::cerr << "Ptable8 const vector psite\n";
    this->reserve(p.size());
    std::for_each( p.cbegin(), p.cend(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, Seq8(p.second,dna_poly_alphabet)) );//poly8::dna2vtype(p.second)) );
		   });
  }
  
  Ptable8::Ptable8( std::vector<polymorphicSite> & p) : base()
  {
    std::cerr << "Ptable8 vector psite\n";
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( polymorphicSite & p ) {
		     this->push_back( poly8site(p.first,Seq8(p.second,std::ref(dna_poly_alphabet))) );//poly8::dna2vtype(p.second)) );
						//poly8::dna2vtype(std::ref(p.second))) );
		   });
    p.clear();
  }
}
