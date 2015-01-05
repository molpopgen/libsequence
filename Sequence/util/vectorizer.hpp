/*!
  \file Sequence/util/vectorizer.hpp
  \brief Definition of Sequence::util::vectorizer
 */

/*! \namespace Sequence::util
  \brief Routines related to low-level C++ stuff
 */

#ifndef __SEQUENCE_UTIL_VECTORIZER_HPP__
#define __SEQUENCE_UTIL_VECTORIZER_HPP__

#include <vector>

namespace Sequence {
  namespace util {
    /*!
      \brief Facilitates safe inheritance from std::vector<T, allocator>
      Privately inherits from std::vector and makes the public interface of that
      class accessible via using declarations.

      Classes may safely inherit publicly from this class with no fear of memory 
      leaks, etc.

      \note The interface is the C++11 interface to std::vector
     */
    template<typename T, typename allocator = std::allocator<T> >
    class vectorizer : private std::vector<T,allocator>
    {
    public:
      using vectorizer_base_t = std::vector<T,allocator>;

      //import base types
      using typename vectorizer_base_t::value_type;
      using typename vectorizer_base_t::allocator_type;
      using typename vectorizer_base_t::reference;
      using typename vectorizer_base_t::const_reference;
      using typename vectorizer_base_t::pointer;
      using typename vectorizer_base_t::const_pointer;
      using typename vectorizer_base_t::iterator;
      using typename vectorizer_base_t::const_iterator;
      using typename vectorizer_base_t::reverse_iterator;
      using typename vectorizer_base_t::const_reverse_iterator;
      using typename vectorizer_base_t::difference_type;
      using typename vectorizer_base_t::size_type;

      //Inherit constructors and assignment operator
      using vectorizer_base_t::vectorizer_base_t;
      using vectorizer_base_t::operator=;
    
      //Iterators:
      using vectorizer_base_t::begin;
      using vectorizer_base_t::end;
      using vectorizer_base_t::rbegin;
      using vectorizer_base_t::rend;
      using vectorizer_base_t::crbegin;
      using vectorizer_base_t::crend;
      using vectorizer_base_t::cbegin;
      using vectorizer_base_t::cend;

      //Capacity:
      using vectorizer_base_t::size;
      using vectorizer_base_t::max_size;
      using vectorizer_base_t::resize;
      using vectorizer_base_t::capacity;
      using vectorizer_base_t::empty;
      using vectorizer_base_t::reserve;
      using vectorizer_base_t::shrink_to_fit;

      //Element access:
      using vectorizer_base_t::operator[];
      using vectorizer_base_t::at;
      using vectorizer_base_t::front;
      using vectorizer_base_t::back;
      using vectorizer_base_t::data;

      //Modifiers:
      using vectorizer_base_t::assign;
      using vectorizer_base_t::push_back;
      using vectorizer_base_t::pop_back;
      using vectorizer_base_t::insert;
      using vectorizer_base_t::erase;
      using vectorizer_base_t::swap;
      using vectorizer_base_t::clear;
      using vectorizer_base_t::emplace_back;
      using vectorizer_base_t::emplace;

      //Allocator:
      using vectorizer_base_t::get_allocator;
    };
  }
}
#endif
    
