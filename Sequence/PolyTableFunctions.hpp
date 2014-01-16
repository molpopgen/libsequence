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

#ifndef __POLYTABLE_FUNCTIONS_HPP__
#define __POLYTABLE_FUNCTIONS_HPP__
#include <Sequence/SeqExceptions.hpp>

/*! \file PolyTableFunctions.hpp
  There are lots of ways one may want to manipulate a Sequence::PolyTable.
  In the past, useful manipulation functions were added as member functions.
  That quickly got out of control, breaking binary compatibility left and 
  right between different library versions.  Declarations for such functions
  are now being added to this file.
  \short Operations on non-const Sequence::PolyTable objects
 */

namespace Sequence
{
  class PolyTable;
  bool containsCharacter(const PolyTable * t,
			 const char & ch);

  void fillIn(PolyTable * t,
	      const unsigned & refseq = 0,
	      const char & identical = '.');

  void addIdentityChar(PolyTable *t,
		       const unsigned & refseq = 0,
		       const char & identical = '.');

  void RemoveGaps(PolyTable *t,
		  const char & gapchar = '-') ;

  void RemoveInvariantColumns(PolyTable *t,
			      const bool & skipOutgroup = false,
			      const unsigned & outgroup = 0) ;

  bool PolyTableValid(const PolyTable * t);
}
#endif
