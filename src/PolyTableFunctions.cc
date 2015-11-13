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

#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/PolyTable.hpp>
#include <algorithm>
#include <set>
#include <cctype>

namespace Sequence
{
  bool containsCharacter(const PolyTable * t,
  			 const char ch)
  {
    for( PolyTable::const_data_iterator itr = t->begin() ;
  	 itr < t->end() ;
  	 ++itr )
      {
  	if ( itr->find(ch) != std::string::npos )
  	  {
  	    return true;
  	  }
      }
    return false;
  }

  bool polyTableValid(const PolyTable * table)
  {
    for ( PolyTable::const_data_iterator itr = table->begin() ;
	  itr < table->end() ; 
	  ++itr )
      {
	if ( (std::find_if(itr->begin(),itr->end(),invalidPolyChar()) != itr->end())
	     || ( itr->length() != table->numsites() ) )
	  {
	    return false;
	  }
      }
    return true;
  }
}
