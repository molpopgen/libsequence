#ifndef __SEQUENCE_SUMMSTATS_UHAPS_HPP__
#define __SEQUENCE_SUMMSTATS_UHAPS_HPP__

#include <Sequence/Seq8.hpp>
#include <Sequence/polySiteVector8.hpp>
#include <Sequence/PolyTable.hpp>
#include <list>
#include <vector>

namespace Sequence 
{
  class Uhaps
  {
  private:
    using ustring_ctr_t = std::list<Seq8>;
    using ustring_const_iterator = ustring_ctr_t::const_iterator;
    using ustring_itr_ctr_t = std::vector< ustring_const_iterator >;
    mutable ustring_ctr_t ustrings;
    ustring_itr_ctr_t ustring_itrs;
    int populate( const polySiteVector8 & );
    int populate( const PolyTable & );
    int populate( PolyTable && );
  public:
    Uhaps() = default;
    Uhaps( const polySiteVector8 & );
    Uhaps( const PolyTable & );
    using iterator = ustring_itr_ctr_t::iterator;
    using const_iterator = ustring_itr_ctr_t::const_iterator;
    using const_reference = ustring_ctr_t::const_reference;
    using size_type = ustring_itr_ctr_t::size_type;
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    const_iterator cbegin() const;
    const_iterator cend() const;
    const_reference operator[]( const size_type & i ) const;
    size_type size() const;
    bool empty() const;
  };
}
#endif
