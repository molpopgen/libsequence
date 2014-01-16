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

/*! \file typedefs.hpp
  Typedefs used in the library are defined here.
  Wherever possible, types from namespace std
  are given forward declarations. 
  @brief typedefs used by libsequence
*/
#ifndef __SEQUENCE_TYPEDEFS_HPP
#define __SEQUENCE_TYPEDEFS_HPP
#include <vector>

namespace Sequence
{
  /*! 
    A CodonUsageTable is a vector of pairs.  In each pair,
    the first element is the codon, and the second element
    is an integer counting the number of occurrences of 
    the codon
  */
  typedef std::vector< std::pair<std::string,int> > CodonUsageTable;

  /*!
    For polymorphism data, a Site can be represented as
    a position (a double) and the characters at 
    that positions (a std::string)
  */
  typedef std::pair< double, std::string > polymorphicSite;

  /*!
    A polymorphism data set can be represented as
    a vector containing a sequence of polymorphicSite
  */
  typedef std::vector< polymorphicSite > polySiteVector;
}
#endif
