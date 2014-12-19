#ifndef __SEQUENCE_PTABLE8_HPP__
#define __SEQUENCE_PTABLE8_HPP__

#include <Sequence/Poly8.hpp>
#include <Sequence/typedefs.hpp>
#include <utility>

namespace Sequence
{
  //! @brief Analog to Sequence::polymorphicSite
  using poly8site = std::pair<double,Sequence::poly8::vtype>;
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
