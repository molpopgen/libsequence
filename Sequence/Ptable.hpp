#ifndef __SEQUENCE_PTABLE_HPP__
#define __SEQUENCE_PTABLE_HPP__

#include <Sequence/typedefs.hpp>

namespace Sequence
{
  //forward declaration
  class PolyTable;

  /*! \class Sequence::Ptable Sequence/Ptable.hpp
    \ingroup seqio
    A "site-major" class to manipulate SNP data.
    \example Ptable_test.cc
   */
  class Ptable : public std::vector<polymorphicSite>
  {
  private:
  public:
    typedef std::vector<polymorphicSite> base;

    //constructors
    Ptable( const std::initializer_list<polymorphicSite> & );
    Ptable( const PolyTable & );
    Ptable( const PolyTable * );
  };
} //ns Sequence

#endif
