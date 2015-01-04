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
  class Ptable8 : private std::vector< poly8site >
  {
  private:
  public:
    using base = std::vector< poly8site >;
    using base::base;
    //Constructors 
    Ptable8() = default;
    Ptable8( const std::initializer_list<polymorphicSite> & );
    Ptable8( const PolyTable & );
    Ptable8( const PolyTable * );
    Ptable8( const std::vector<polymorphicSite> & );
    Ptable8( std::vector<polymorphicSite> & );

    //import base types
    //import base functions
    using base::size;
    using base::max_size;
    using base::capacity;
    using base::erase;
    using base::clear;
    using base::insert;
    using base::swap;
    using base::pop_back;
    using base::reserve;
    using base::empty;
    using base::resize;
    using base::begin;
    using base::end;
    using base::rbegin;
    using base::rend;
    using base::crbegin;
    using base::crend;
    using base::cbegin;
    using base::cend;
    using base::operator[];
    using base::operator=;
    using base::at;
    using base::push_back;
    using base::emplace_back;
    using base::front;
    using base::back;
    using base::data;
  };
}

#endif
