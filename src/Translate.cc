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

#include <Sequence/Translate.hpp>
#include <algorithm>
#include <cctype>

namespace
{
  char Universal (const char codon[4],
		  const char & gapchar)
  {
    //handle gaps A codon with a single or 2
    //gap characters returns an X, since it's ambiguous,
    //and a codon that is all gaps (---)
    //returns a -
    std::ptrdiff_t ngaps = (std::count(&codon[0],&codon[0]+3,gapchar));
    if (ngaps == 3)
      {
	return '-';
      }
    else if (ngaps > 0 && ngaps < 3)
      {
	return 'X';
      }

    switch (codon[0])
      {
        //first position is A
      case 'A':
        switch (codon[1])
          {
          case 'T':
            switch (codon[2])
              {
              case 'G':
                return 'M';
                break;
              default:
                return 'I';
                break;
              }
            break;
          case 'C':
            return 'T';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'T':
                return 'N';
                break;
              case 'C':
                return 'N';
                break;
              case 'A':
                return 'K';
                break;
              case 'G':
                return 'K';
                break;
              }
          case 'G':
            switch (codon[2])
              {
              case 'A':
                return 'R';
                break;
              case 'G':
                return 'R';
                break;
              default:
                return 'S';
                break;
              }
            break;
          }
        break;
        //first Position is T
      case 'T':
        switch (codon[1])
          {
          case 'T':
            switch (codon[2])
              {
              case 'T':
                return 'F';
                break;
              case 'C':
                return 'F';
                break;
              case 'A':
                return 'L';
                break;
              case 'G':
                return 'L';
                break;
              }
            break;
          case 'C':
            return 'S';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'T':
                return 'Y';
                break;
              case 'C':
                return 'Y';
                break;
              default:
                return '*';
                break;
              }
            break;
          case 'G':
            switch (codon[2])
              {
              case 'A':
                return '*';
                break;
              case 'G':
                return 'W';
                break;
              default:
                return 'C';
                break;
              }
            break;
          }

        break;
        //first is G
      case 'G':
        switch (codon[1])
          {
          case 'T':
            return 'V';
            break;
          case 'C':
            return 'A';
            break;
          case 'G':
            return 'G';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'A':
                return 'E';
                break;
              case 'G':
                return 'E';
                break;
              default:
                return 'D';
                break;
              }
            break;
          }
        break;
        //first is C
      case 'C':
        switch (codon[1])
          {
          case 'T':
            return 'L';
            break;
          case 'C':
            return 'P';
            break;
          case 'G':
            return 'R';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'A':
                return 'Q';
                break;
              case 'G':
                return 'Q';
                break;
              default:
                return 'H';
                break;
              }
            break;
          }
        break;
      default:
        return 'X';
        break;
      }
    return 'X';
  }
}

namespace Sequence
{
  std::string Translate(std::string::const_iterator beg,
			std::string::const_iterator end,
			Sequence::GeneticCodes genetic_code,
			const char & gapchar)
    

  {
    if (beg > (end-3)) //if the range is less than 3 in length (1 codon), return an empty string
      return std::string();

    std::string translation;
    char codon[4];
    codon[3] = '\0'; //null terminate for cleanliness
    while(beg < end)
      {
        codon[0] = char(std::toupper(*beg));
        codon[1] = char(std::toupper(*(beg+1)));
        codon[2] = char(std::toupper(*(beg+2)));
        switch (genetic_code)
          {
          case Sequence::UNIVERSAL:
            translation += Universal (codon,gapchar);
            break;
          default:
            throw SeqException ("Translate.cc: Translate(), invalid genetic code passed");
            break;
          }
        beg += 3;
      }
    return translation;
  }

}
