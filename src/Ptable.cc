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
}
