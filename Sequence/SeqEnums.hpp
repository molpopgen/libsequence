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

#ifndef SEQENUMS_H
#define SEQENUMS_H
#include <cstdint>
/*! \file SeqEnums.hpp
  Defines a handfull of enumeration types useful
  for sequence data.
  @brief Definition of enumeration types
*/

namespace Sequence
  {
    /*! \enum Sequence::GeneticCodes
      Only UNIVERSAL (= 0)  is currently supported.
      The order of the genetic codes is that of NCBI's code tables, available at 
      http://www.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=c#SG2\n
    */
    enum class GeneticCodes : std::int16_t {UNIVERSAL};
    /*! \enum Sequence::Mutations
      Values: Unknown=0,Ts, and Tv.\n
      Unknown means unknown, Ts means transition, Tv means transversion
    */
    enum class Mutations : std::int8_t {Unknown,Ts,Tv};
}
#endif
