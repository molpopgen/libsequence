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

//  -*- C++ -*- 
#ifndef __COUNTING_OPERATORS_TCC__
#define __COUNTING_OPERATORS_TCC__
#include <Sequence/CountingOperators.hpp>

namespace Sequence
{
  template<typename key, typename value>
  std::vector<std::pair<key,value> >
  operator+(const std::vector<std::pair<key,value> > &lhs,
	    const std::vector<std::pair<key,value> > &rhs)
  {
    typedef std::vector<std::pair<key,value> > rt; //return type
    rt rv(lhs);
    for(typename rt::const_iterator itr = rhs.begin();
	itr < rhs.end();
	++itr)
      {
	typename rt::iterator i = std::find_if(rv.begin(),
					       rv.end(),
					       std::bind2nd(first_is_equal<key,value>(),*itr));
	if (i != rv.end())
	  {
	    i->second += itr->second;
	  }
	else
	  {
	    rv.push_back(*itr);
  }
      }
    return rv;
  }

  template<typename key, typename value>
  std::vector<std::pair<key,value> >
  operator+=( std::vector<std::pair<key,value> > &lhs,
	      const std::vector<std::pair<key,value> > &rhs)
  {
    return lhs=lhs+rhs;
  }

  template< typename key, typename value,typename comparison>
  std::map<key,value,comparison> operator+(const std::map<key,value,comparison> &lhs,
			  const std::map<key,value,comparison> &rhs)
  {
    typedef std::map<key,value,comparison> rt;
    rt rv(lhs);
    for(typename rt::const_iterator itr = rhs.begin() ;
	itr != rhs.end() ;
	++itr)
      {
	typename rt::iterator i = rv.find(itr->first);
	if ( i != rv.end() )
	  {
	    i->second += itr->second;
	  }
	else
	  {
	    rv[itr->first] = itr->second;
	  }
      }
    return rv;
  }

  template< typename key, typename value,typename comparison>
  std::map<key,value,comparison> operator+=( std::map<key,value,comparison> &lhs,
			    const std::map<key,value,comparison> &rhs)
  {
    return lhs=lhs+rhs;
  }
}
#endif
