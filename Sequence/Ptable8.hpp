#ifndef __SEQUENCE_PTABLE8_HPP__
#define __SEQUENCE_PTABLE8_HPP__

#include <Sequence/Seq8.hpp>
#include <Sequence/polySiteVector.hpp>
#include <utility>

namespace Sequence
{
  //! @brief Analog to Sequence::polymorphicSite
  using poly8site = std::pair<double,Seq8>;
  //fwd declarations
  class PolyTable;
  class Ptable8 : public std::vector< poly8site >
  {
  private:
  public:
    using base = std::vector< poly8site >;

    //Constructors 
    Ptable8() = default;
    Ptable8( const std::initializer_list<polymorphicSite> & );
    Ptable8( const PolyTable & );
    Ptable8( const PolyTable * );
    Ptable8( const std::vector<polymorphicSite> & );
    Ptable8( std::vector<polymorphicSite> & );
  };
}

#endif
