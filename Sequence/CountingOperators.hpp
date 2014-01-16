/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

*/

/*!
  \defgroup operators Operator overloads
  \short Operator overloads
*/
/*!
  \file CountingOperators.hpp
  \brief Declarations of operators to add associative containers together

  A lot of biological computation involves counting.  Associative containers
  such as vectors of pairs or maps (or hashes if you've got them) are natural
  data structures to keeps track of counts.  However, keeping track of a running
  total of counts (say codon usage summed accross different sequences) is hampered
  by the lack of addition operators for associative containers.  This file 
  defines template operator+ and += for std::vector<std::pair<key,I> > and 
  std::map<key,I>.  The results of these operations are best explained by example:
  \code
  #include <Sequence/CountingOperators.hpp>
  #include <map>
  std::map<char,unsigned> baseCounts,baseCounts2;
  baseCounts['A'] = 5;
  baseCounts['G'] = 10;
  baseCounts2['A'] = 11;
  baseCounts2['C'] = 17;
  std::map<char,unsigned> baseCounts3 = baseCounts + baseCounts2;
  for (std::map<char,unsigned>::const_iterator i = baseCounts3.begin() ;
  i != baseCounts3.end() ;
  ++i)
  {
  std::cout << i->first << '\t' << i->second << '\n';
  }
  \endcode
  The output of the above will be (possibly in a different order):
  A  16\n
  C  17\n
  G  10\n

  Therefore, the result of the addition operation is an associative container
  containing all the "key" elements of the 2 containers being added.  And when
  those containers contain common keys, the values associated with those keys
  are summed.

  These operators are not restricted in terms of types for the template arguments.
  The only requirements are that operator== is defined for the key, and operator+
  must be valid for the value type.
*/
#ifndef __COUNTING_OPERATORS_HPP__
#define __COUNTING_OPERATORS_HPP__
#include <algorithm>
#include <functional>
#include <map>
#include <set>
namespace Sequence
{
  template<typename key,typename value>
  struct first_is_equal : public std::binary_function< std::pair<key,value>,
						       std::pair<key,value>,bool>
			  /*! \struct first_is_equal Sequence/CountingOperators.hpp
			    \brief Functor that checks for equality of first member of two pairs
			   */
  {
    inline bool operator()(const std::pair<key,value> &l,
			   const std::pair<key,value> &r) const
    /*!
      \return true if l.first == r.first, false otherwise
    */
    {
      return l.first==r.first;
    }
  };
  
  /*!
    Add 2 vectors of pairs together.
    The return vector contains all the elements of type T
    present in \a lhs and \a rhs.  For all elements T
    that \a lhs and \a rhs have in common, the associated
    value of element I is the sum of the values in \a lhs
    and \a rhs.  
    \note operator== must be defined for type T, and operator+ must
    be valid for type value
    \ingroup operators
  */
  template<typename key, typename value>
  std::vector<std::pair<key,value> >
  operator+(const std::vector<std::pair<key,value> > &lhs,
	    const std::vector<std::pair<key,value> > &rhs);


  /*!
    operator+= for two vectors of pairs.
    \note operator== must be defined for type key, and operator+ must
    be valid for type value
    \ingroup operators
  */
  template<typename key, typename value>
  std::vector<std::pair<key,value> >
  operator+=( std::vector<std::pair<key,value> > &lhs,
	      const std::vector<std::pair<key,value> > &rhs);

  /*!
    Add 2 maps together.
    The return vector contains all the elements of type T
    present in \a lhs and \a rhs.  For all elements T
    that \a lhs and \a rhs have in common, the associated
    value of element I is the sum of the values in \a lhs
    and \a rhs.
    \note operator== must be defined for type T, and operator+ must
    be valid for type value
    \ingroup operators
  */
  template< typename key, typename value,
	    typename comparison>
  std::map<key,value,comparison> operator+(const std::map<key,value,comparison> &lhs,
					   const std::map<key,value,comparison> &rhs);

  /*!
    operator+= for two maps
    \note operator== must be defined for type key, and operator+ must
    be valid for type value
    \ingroup operators
  */
  template< typename key, typename value, typename comparison>
  std::map<key,value,comparison> operator+=( std::map<key,value,comparison> &lhs,
			    const std::map<key,value,comparison> &rhs);
}
#include <Sequence/bits/CountingOperators.tcc>
#endif
