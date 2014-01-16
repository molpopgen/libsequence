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

#include <Sequence/ComplementBase.hpp>
#include <cctype>


namespace Sequence
  {
  void ComplementBase::operator()(char &ch) const
  /*!
    Converts ch to the complement of ch. Use with std::for_each to complement a range
   */
    {
      //      char uch = std::toupper(ch);
      switch(ch)
        {
        case 'A':
          ch = 'T';
          break;
	case 'a':
	  ch = 't';
	  break;
        case 'T':
          ch = 'A';
          break;
	case 't':
	  ch = 'a';
	  break;
        case 'G':
          ch = 'C';
          break;
	case 'g':
	  ch = 'c';
	  break;
        case 'C':
          ch = 'G';
          break;
	case 'c':
	  ch = 'g';
	  break;
        case 'M':
          ch = 'K';
          break;
	case 'm':
	  ch ='k';
	  break;
        case 'K':
          ch = 'M';
          break;
	case 'k':
	  ch = 'm';
        case 'R':
          ch = 'Y';
          break;
	case 'r':
	  ch = 'y';
	  break;
        case 'Y':
          ch = 'R';
          break;
	case 'y':
	  ch = 'r';
	  break;
        case 'W':
          ch = 'W';
          break;
	case 'w':
	  ch = 'w';
	  break;
        case 'S':
          ch = 'S';
          break;
	case 's':
	  ch = 's';
	  break;
        case 'B':
          ch = 'V';
          break;
	case 'b':
	  ch = 'v';
	  break;
        case 'V':
          ch = 'B';
          break;
	case 'v':
	  ch = 'b';
	  break;
        case 'D':
          ch = 'H';
          break;
	case 'd':
	  ch = 'h';
	  break;
        case 'H':
          ch = 'D';
          break;
	case 'h':
	  ch = 'd';
	  break;
        case 'X':
          ch = 'X';
          break;
	case 'x':
	  ch = 'x';
	  break;
        case 'N':
          ch = 'N';
          break;
	case 'n':
	  ch = 'n';
        case '-':
          ch = '-';
          break;
        default:
          ch = '?';
          break;
        }
    }

}

