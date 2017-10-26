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

#ifndef __POLYSITEVECTOR_MANIP_HPP__
#define __POLYSITEVECTOR_MANIP_HPP__


/*! \file Sequence/polySiteVector.hpp
  @brief Site-major variation tables in ASCII format
*/

#include <string>
#include <utility>
#include <vector>

namespace Sequence
{
  class PolyTable;

  /*!
    For polymorphism data, a Site can be represented as
    a position (a double) and the characters at 
    that positions (a std::string)
  */
  using polymorphicSite = std::pair< double, std::string >__attribute__((deprecated));

  /*!
    A polymorphism data set can be represented as
    a vector containing a sequence of polymorphicSite
  */
  using polySiteVector = std::vector< polymorphicSite >__attribute__((deprecated));

  polySiteVector make_polySiteVector(const Sequence::PolyTable & data)__attribute__((deprecated));
}
#endif
