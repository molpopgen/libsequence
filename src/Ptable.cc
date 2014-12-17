#include <Sequence/Ptable.hpp>
#include <Sequence/PolyTable.hpp>

namespace Sequence
{
  Ptable::Ptable( const PolyTable & pt ) : base(pt.sbegin(),pt.send())
  {
  }
  
  Ptable::Ptable( const PolyTable * pt) : base(pt->sbegin(),pt->send())
  {
  }

  Ptable::Ptable( const std::initializer_list<polymorphicSite> & il) : base(il.begin(),il.end())
  {
  }

  Ptable::Ptable( const std::vector<polymorphicSite> & v ) : base(v)
  {
  }

  Ptable::Ptable( std::vector<polymorphicSite> && v ) : base(std::move(v))
  {
  }
}
