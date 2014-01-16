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

#include <string>
#include <functional>
#include <cctype>

namespace Sequence
{
  struct lt_nocase : public std::binary_function<char, char, bool> 
  /*
    Fast, clean, and slightly incorrect implementations of case-insensitive
    character comparisons
    \warning not locale-aware!
   */
  { 
    bool operator()(char x, char y) const 
    { 
      return std::toupper(static_cast<unsigned char>(x)) 
	< std::toupper(static_cast<unsigned char>(y)); 
    } 
  }; 
}
