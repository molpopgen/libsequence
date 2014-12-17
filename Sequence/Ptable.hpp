//! \file Sequence/Ptable.hpp @brief Site-major variation tables
#ifndef __SEQUENCE_PTABLE_HPP__
#define __SEQUENCE_PTABLE_HPP__

#include <Sequence/typedefs.hpp>

namespace Sequence
{
  //forward declaration
  class PolyTable;

  /*! \class Sequence::Ptable Sequence/Ptable.hpp
    \ingroup polytables
    A "site-major" class to manipulate SNP data.
    \example Ptable_test.cc
   */
  class Ptable : public std::vector<polymorphicSite>
  {
  private:
  public:
    using base = std::vector<polymorphicSite>;

    //constructors
    //! Empty
    Ptable() : std::vector<polymorphicSite>() {}
    //! Construct from an initializer_list
    Ptable( const std::initializer_list<polymorphicSite> & );
    //! Construct from const PolyTable
    Ptable( const PolyTable & );
    //! Construct from const PolyTable *
    Ptable( const PolyTable * );
    Ptable( const std::vector<polymorphicSite> & );
    Ptable( std::vector<polymorphicSite> && );
  };
} //ns Sequence

#endif
