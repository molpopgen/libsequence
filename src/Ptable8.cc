#include <Sequence/Ptable8.hpp>
#include <Sequence/PolyTable.hpp>
namespace Sequence
{
  Ptable8::Ptable8( const std::initializer_list<polymorphicSite> & p) : base()
  {
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, poly8::dna2vtype(p.second)) );
		   });
  }
  
  Ptable8::Ptable8( const PolyTable & p) : base()
  {
    this->reserve(p.numsites());
    std::for_each( p.sbegin(), p.send(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, poly8::dna2vtype(p.second)) );
		   });
  }

  Ptable8::Ptable8( const PolyTable * p) : base()
  {
    *this = Ptable8(*p);
  }

  Ptable8::Ptable8( const std::vector<polymorphicSite> & p) : base()
  {
    this->reserve(p.size());
    std::for_each( p.cbegin(), p.cend(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, poly8::dna2vtype(p.second)) );
		   });
  }
  
  Ptable8::Ptable8( std::vector<polymorphicSite> & p) : base()
  {
    this->reserve(p.size());
    std::for_each( p.begin(), p.end(),
		   [this]( const polymorphicSite & p ) {
		     this->push_back( poly8site(p.first, poly8::dna2vtype(std::ref(p.second))) );
		   });
    p.clear();
  }
}
