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
#include <limits>

namespace Sequence
{
  WeightingScheme2::weights2_t Unweighted2::operator()(const std::string &, const std::string &,Sequence::GeneticCodes) const
  /*!
    Calculate actually calculates the weights for each branch
    \param codon1 a std::string of length 3 representing a sense codon
    \param codon2 a std::string of length 3 representing a sense codon
  */
  {
    return weights2_t({1.,1.});
  }

  WeightingScheme3::weights3_t Unweighted3::operator()(const std::string &, const std::string &,Sequence::GeneticCodes ) const
  /*!
    Calculate actually calculates the weights for each branch
    \param codon1 a std::string of length 3 representing a sense codon
    \param codon2 a std::string of length 3 representing a sense codon
  */
  {
    return weights3_t({1.,1.,1.,1.,1.,1.,});
  }
}
