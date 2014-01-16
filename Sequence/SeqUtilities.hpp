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

#ifndef __SEQ_UTILITIES_H__
#define __SEQ_UTILITIES_H__
#include <iterator>
#include <algorithm>
#include <map>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>

/*! \file SeqUtilities.hpp
  This file declares various functions (mostly templates) that may be useful in programs.
  @short Declaration of Sequence::makeCountList (an alternative to Sequence::stateCounter),
  Sequence::internalGapCheck
*/

/*! \fn  std::map< typename std::iterator_traits< Iterator >::value_type, unsigned > Sequence::makeCountList(Iterator beg, Iterator end) 
  \param beg an iterator
  \param end an iterator
  \return a std::map< type, unsigned >, where type is the
  iterator_traits<Iterator>::value_type of Iterator.  
  The keys are the (unique) elements
  present in the range, and the unsinged values the numbers
  of times each element occurs occur
  \note This function can be used as an alternative to
  Sequence::stateCounter if you want to count more than just strict DNA
  characters.  
*/

/*! \fn  bool Sequence::internalGapCheck(Iterator beg, Iterator end,const char & gapchar = '-',const unsigned & mod = 3)
  This function checks a range for internal gaps that meet a certain
  length requirement.  The requirement is that length%mod == 0.  The value 
  true is returned if this is not the case, false otherwise.  One use of
  this function may be to check that the internal gaps in an aligned cds 
 sequence are all multiples of 3 in length.
*/
namespace Sequence
{
  template<typename Iterator> 
  std::map<typename std::iterator_traits<Iterator>::value_type,unsigned>
  makeCountList( Iterator beg,  Iterator end)
  {
    typedef typename std::iterator_traits<Iterator>::value_type type;
    typedef typename std::map<type,unsigned> maptype;
    typedef typename maptype::iterator mitr;
    maptype m;
    mitr i;
    while(beg < end)
      {
	i = m.find(*beg);
	if ( i != m.end() )
	  {
	    (*i).second++;
	  }
	else
	  {
	    m[*beg] = 1;
	  }
	beg++;
      }
    return m;
  }

  template<typename Iterator>
  bool
  internalGapCheck(Iterator beg, Iterator end,
		   const char & gapchar = '-',
		   const unsigned & mod = 3)
  {
    //the value type of the iterator must be convertible to a char
    BOOST_STATIC_ASSERT( (boost::is_convertible<typename std::iterator_traits<Iterator>::value_type,char>::value) );
    if(beg==end) return false;
    Iterator itr,effective_end=end-1;
    //find the first non-gap at the end of range
    while(effective_end > beg)
      {
	if (*effective_end != gapchar)
	  {
	    break;
	  }
	else
	  {
	    --effective_end;
	  }
      }
    itr = std::find_if(beg,effective_end,
		       std::bind2nd(std::not_equal_to<char>(),gapchar));
    ++itr;
    unsigned cont_gap = 0;
    while ( itr < effective_end )
      {
	if ( *itr == gapchar )
	  {
	    cont_gap++;
	  }
	else
	  {
	    if (cont_gap % mod != 0.)
	      {
		//gap fails our length requirement
		return true;
	      }
	    else
	      {
		cont_gap = 0;
	      }
	  }
	itr++;
      }
    return false;
  }
}//namespace Sequence

#endif
