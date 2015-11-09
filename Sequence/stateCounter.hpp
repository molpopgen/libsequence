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

#ifndef _STATE_COUNTER_H_
#define _STATE_COUNTER_H_
#include <functional>
/*! \file stateCounter.hpp
@brief declaration of Sequence::stateCounter, a class to keep track of nucleotide counts either at a site
in an alignment, or along a sequence
*/
/*!
  \class Sequence::stateCounter Sequence/stateCounter.hpp
  \ingroup functors
  \warning class data are public.  Use responsibly.
  @short keep track of state counts at a site in an alignment or along a sequence
*/
namespace Sequence
  {
    class stateCounter : public std::unary_function<char,void>
    {
    private:
      char _gap;
    public:
      typedef unsigned size_type;
      size_type a,g,c,t,zero,one,gap,n;
      bool ndna;
      stateCounter(char gapchar = '-');
      stateCounter(stateCounter &&) = default;
      stateCounter(const stateCounter &) = default;
      stateCounter & operator=(const stateCounter &)=default;
      stateCounter & operator=( stateCounter &&)=default;
      void operator()(const char &ch);
      size_type nStates(void) const;
    };
}
#endif
