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

#include <Sequence/Unweighted.hpp>

namespace Sequence
  {
  void Unweighted2::Calculate(const std::string &codon1, const std::string &codon2) const
  /*!
    Calculate actually calculates the weights for each branch
    \param codon1 a std::string of length 3 representing a sense codon
    \param codon2 a std::string of length 3 representing a sense codon
  */
    {
      __weights[0] = 0.5;
      __weights[1] = 0.5;
    }

  double* Unweighted2::weights(void) const
  /*!
    \return a double * of size 2 (1 value for each branch)
  */
    {
      return __weights;
    }

  void Unweighted3::Calculate(const std::string &codon1, const std::string &codon2) const
  /*!
    Calculate actually calculates the weights for each branch
    \param codon1 a std::string of length 3 representing a sense codon
    \param codon2 a std::string of length 3 representing a sense codon
  */
    {
      __weights[0] = 1.0/6.0;
      __weights[1] = 1.0/6.0;
      __weights[2] = 1.0/6.0;
      __weights[3] = 1.0/6.0;
      __weights[4] = 1.0/6.0;
      __weights[5] = 1.0/6.0;
    }

  double* Unweighted3::weights(void) const
  /*!
    \return a double * of size 6 (1 value for each branch)
  */
    {
      return __weights;
    }
}
