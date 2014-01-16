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

#ifndef __PATHWAYHELPER_H__
#define __PATHWAYHELPER_H__
/*! \file PathwayHelper.hpp
  @brief declarations of Sequence::Intermediates2 and Sequence::Intermediates3
*/
/*!
  \defgroup CodonPaths Classes and functions to aid in the calculations of the pathways between two codons
  This group of classes and functions deals with determining 
  either the counts of silent and replacement differences between codons
  or the intermedate codons that occurs between two different codons
*/
#include <string>

namespace Sequence
{
  void Intermediates2(std::string *intermediates,const std::string &codon1, const std::string &codon2);
  void Intermediates3(std::string *intermediates, const  std::string &codon1, const std::string &codon2);
}
#endif

