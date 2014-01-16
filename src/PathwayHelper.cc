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
#include <vector>
#include <cassert>
#include <cctype>

using std::string;

namespace Sequence
{
  void Intermediates2(string *intermediates,
		      const std::string &codon1,
		      const std::string &codon2)
  /*!
    \param intermediates a string[2] in which we will place the intermediate codons
    \param codon1 a codon
    \param codon2 a codon
    \short Calculate the intermediate codons between a pair of codons diverged at 2 positions
    \ingroup CodonPaths
  */
  {
    intermediates[0].resize(3);
    intermediates[1].resize(3);

    unsigned i,j;
    int pos[2];

    for (i = 0, j = 0; i <= 2; ++i)
      if (char(std::toupper(codon1[i])) != char(std::toupper(codon2[i])))
	pos[j++] = i;
    switch (pos[0])
      {
      case 0:
	intermediates[0][0] = char(std::toupper(codon2[0]));
	intermediates[0][1] = char(std::toupper(codon1[1]));
	intermediates[0][2] = char(std::toupper(codon1[2]));
	break;
      case 1:
	intermediates[0][0] = char(std::toupper(codon1[0]));
	intermediates[0][1] = char(std::toupper(codon2[1]));
	intermediates[0][2] = char(std::toupper(codon1[2]));
	break;
      case 2:
	intermediates[0][0] = char(std::toupper(codon1[0]));
	intermediates[0][1] = char(std::toupper(codon1[1]));
	intermediates[0][2] = char(std::toupper(codon2[2]));
	break;
      }

    switch(pos[1])
      {
      case 0:
	intermediates[1][0] = char(std::toupper(codon2[0]));
	intermediates[1][1] = char(std::toupper(codon1[1]));
	intermediates[1][2] = char(std::toupper(codon1[2]));
	break;
      case 1:
	intermediates[1][0] = char(std::toupper(codon1[0]));
	intermediates[1][1] = char(std::toupper(codon2[1]));
	intermediates[1][2] = char(std::toupper(codon1[2]));
	break;
      case 2:
	intermediates[1][0] = char(std::toupper(codon1[0]));
	intermediates[1][1] = char(std::toupper(codon1[1]));
	intermediates[1][2] = char(std::toupper(codon2[2]));
	break;
      }
  }

  void Intermediates3(string *intermediates,const std::string &codon1, const std::string &codon2)
  /*!
    \param intermediates a string[9] in which we will place the intermediate codons
    \param codon1 a codon
    \param codon2 a codon
    \note the storage of the intermediate codons follows the illustration in the documentation of Sequence::ThreeSubs
    \short Calculate the intermediate codons between a pair of codons diverged at 3 positions
    \ingroup CodonPaths
  */
  {
    for(int i = 0 ; i < 9 ;++i)
      intermediates[i].resize(3);

    intermediates[0][0] = char(std::toupper(codon2[0]));
    intermediates[0][1] = char(std::toupper(codon1[1]));
    intermediates[0][2] = char(std::toupper(codon1[2]));

    intermediates[1][0] = char(std::toupper(codon2[0]));
    intermediates[1][1] = char(std::toupper(codon2[1]));
    intermediates[1][2] = char(std::toupper(codon1[2]));

    intermediates[2][0] = char(std::toupper(codon2[0]));
    intermediates[2][1] = char(std::toupper(codon1[1]));
    intermediates[2][2] = char(std::toupper(codon2[2]));

    intermediates[3][0] = char(std::toupper(codon1[0]));
    intermediates[3][1] = char(std::toupper(codon2[1]));
    intermediates[3][2] = char(std::toupper(codon1[2]));

    intermediates[4][0] = char(std::toupper(codon2[0]));
    intermediates[4][1] = char(std::toupper(codon2[1]));
    intermediates[4][2] = char(std::toupper(codon1[2]));

    intermediates[5][0] = char(std::toupper(codon1[0]));
    intermediates[5][1] = char(std::toupper(codon2[1]));
    intermediates[5][2] = char(std::toupper(codon2[2]));

    intermediates[6][0] = char(std::toupper(codon1[0]));
    intermediates[6][1] = char(std::toupper(codon1[1]));
    intermediates[6][2] = char(std::toupper(codon2[2]));

    intermediates[7][0] = char(std::toupper(codon2[0]));
    intermediates[7][1] = char(std::toupper(codon1[1]));
    intermediates[7][2] = char(std::toupper(codon2[2]));

    intermediates[8][0] = char(std::toupper(codon1[0]));
    intermediates[8][1] = char(std::toupper(codon2[1]));
    intermediates[8][2] = char(std::toupper(codon2[2]));
  }
}
