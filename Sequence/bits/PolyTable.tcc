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

// Code for the -*- C++ -*- namespace Sequence::PolyTable template members
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
    BOOST_STATIC_ASSERT( (boost::is_convertible<numeric_type,double>::value) );
    //The character type must be eithe char * or std::string
    BOOST_STATIC_ASSERT( (boost::is_same<string_type,char *>::value  ||
			  boost::is_same<string_type,std::string>::value ) );
  
    positions.resize(_num_positions);
    data.resize(_num_individuals);
    positions.assign(_positions,_positions+_num_positions);
    data.assign(_data,_data+_num_individuals);
    non_const_access = true;
    for(std::vector<std::string>::const_iterator itr = data.begin() ;
	itr < data.end() ; ++itr)
      {
	if (itr->length() != _num_positions)
	  {
	    positions.clear();
	    data.clear();
	    return false;
	  }
      }
    return true;
  }

//   template<typename iterator>
//   bool PolyTable::rear_insert( const iterator beg,
// 			       const iterator end )
//   /*!
//     Insert a range of string types (chr * or std::string)
//     at the end of a PolyTable object.
//    */
//   {
//     typedef typename std::iterator_traits<iterator>::value_type vtype;
//     BOOST_STATIC_ASSERT( (boost::is_same<vtype,char *>::value  ||
// 			  boost::is_same<vtype,std::string>::value ) );
//     data.insert(data.end(),beg,end);
//     return true;
//   }

}

#endif
