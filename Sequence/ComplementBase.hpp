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

#ifndef __COMPLEMENT_BASE_H__
#define __COMPLEMENT_BASE_H__
#include <functional>

/*! \file ComplementBase.hpp
  @brief Delcaration of Sequence::ComplementBase, a function object to return the complement of a DNA nucleotide
*/
/*! 
  \struct Sequence::ComplementBase Sequence/ComplementBase.hpp
  \ingroup functors
  a functor to complement a sequence\n
  example use:
  \code
  //reverse and complement a std::string
  #include <string>
  #include <algorithm>
  #include <Sequence/SeqFunctors.hpp>
  
  int main ()
  {
  std::string seq;
  //fill seq with DNA characters
  std::reverse(seq.begin(),seq.end());
  std::for_each(seq.begin(),seq.end(),Sequence::ComplementBase());
  }
  \endcode
*/
namespace Sequence
  {
  struct ComplementBase : public std::unary_function<char,void>
  {
    void operator()(char &ch) const;
  };
}
#endif
