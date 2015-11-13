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

#include <Sequence/PolyTable.hpp>
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
  //! \return true if \a t contains \a ch, false otherwise
  bool containsCharacter(const PolyTable * t,const char ch);

  /*!
    \return true if the following conditions are met : First, the length of every
    row in the table is equal to the length of the vector of positions.  Second,
    all of the characters in the data are members of the set {A,G,C,T,N,-} 
    (case-insensitive).

    This function is useful if you play around with PolyTable objects in
    non-const contexts, or read them in from files and need to check that
    the data are compatible with other routines in this library.  This 
    routine can be thought of as a PolyTable equivalent to Alignment::validForPolyAnalysis,
    which works on ranges of Sequence::Seq objects.
  */
  bool polyTableValid(const PolyTable * t);

  template<typename T> T copyPolyTable(const T & t)
  {
    return T(std::vector<double>(t.pbegin(),t.pend()),
	     std::vector<std::string>(t.begin(),t.end()));
  }
  
  template<typename T,typename F> T removeColumns( const T & t, const F & f, const bool skipAnc = false, const unsigned anc = 0,const char gapchar = '-' );
  template<typename T> T removeGaps( const T & t, const bool skipAnc = false, const unsigned anc = 0,const char gapchar = '-' );
  template<typename T> T removeInvariantPos(const T & t, const bool skipAnc = false, const unsigned anc = 0,
					    const char gapchar = '-');
  template<typename T> T removeAmbiguous(const T & t, const bool skipAnc = false, const unsigned anc = 0,
					 const char gapchar = '-');
  template<typename T> T removeMissing(const T & t, const bool skipAnc = false, const unsigned anc = 0,
				       const char gapchar = '-');
  template<typename T> T removeMultiHits(const T & t, const bool skipAnc = false, const unsigned anc = 0,
					 const char gapchar = '-');
  template<typename T> T polyTableToBinary(const T & t, const unsigned ref = 0, const char gapchar = '-');
  template<typename T> T polyTableFreqFilter(const T & t, const unsigned mincount,const bool skipAnc = false, const unsigned anc = 0 , const char gapchar = '-');
}
#include <Sequence/bits/PolyTableFunctions.tcc>
#endif
