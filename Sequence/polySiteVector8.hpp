#ifndef __SEQUENCE_POLYSITEVECTOR8_HPP__
#define __SEQUENCE_POLYSITEVECTOR8_HPP__

#include <Sequence/util/vectorizer.hpp>
#include <Sequence/Seq8.hpp>
#include <Sequence/polySiteVector.hpp>
#include <utility>

namespace Sequence
{
  //! @brief Analog to Sequence::polymorphicSite
  using polymorphicSite8 = std::pair<double,Seq8>;
  //fwd declarations
  class PolyTable;
  class polySiteVector8 final : public util::vectorizer<polymorphicSite8> 
  {
  private:
  public:
    using base = util::vectorizer<polymorphicSite8>;
    using base::base;
    //Constructors 
    polySiteVector8() = default;
    polySiteVector8( const std::initializer_list<polymorphicSite> & );
    polySiteVector8( const PolyTable & );
    polySiteVector8( const polySiteVector & );
    polySiteVector8( polySiteVector & );
  };
}

#endif
