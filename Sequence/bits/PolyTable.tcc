// Code for the -*- C++ -*- namespace Sequence::PolyTable template members

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


#ifndef __POLY_TABLE_TCC__
#define __POLY_TABLE_TCC__

#include <iterator>

namespace Sequence
{
  template<typename numeric_type,
	   typename string_type>
  bool PolyTable::assign( const numeric_type * _positions, 
			  const size_t & _num_positions,
			  const string_type * _data,
			  const size_t & _num_individuals )
  { 
    //The numeric array must be convertible to double
    static_assert( std::is_convertible<numeric_type,double>::value,
			"numeric_type must be convertible to double");
    //The character type must be eithe char * or std::string
    static_assert( (std::is_same<string_type,char*>::value || 
		    std::is_same<string_type,std::string>::value),
		   "string_type must be char * or std::string");
  
    first.resize(_num_positions);
    second.resize(_num_individuals);
    first.assign(_positions,_positions+_num_positions);
    second.assign(_data,_data+_num_individuals);
    non_const_access = true;
    for(std::vector<std::string>::const_iterator itr = second.begin() ;
	itr < second.end() ; ++itr)
      {
	if (itr->length() != _num_positions)
	  {
	    first.clear();
	    second.clear();
	    return false;
	  }
      }
    return true;
  }
}

#endif
