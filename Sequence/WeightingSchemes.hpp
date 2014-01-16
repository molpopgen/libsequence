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

#ifndef __WEIGHTING_SCHEMES_H__
#define __WEIGHTING_SCHEMES_H__

#include <string>

/*! \file WeightingSchemes.hpp
  \short abstract interface to weighting schemes when codons differ at 2 positions
 */
/*!
  \defgroup weights Classes and functions relating to pathway weighting for Ka/Ks
  \ingroup kaks
*/
/*!
  \class Sequence::WeightingScheme2 Sequence/WeightingSchemes.hpp
  \ingroup weights
  \warning Before deriving classes from this base, please study the implementations
  of Sequence::TwoSubs.  You should also study Sequence::GranthamWeights2 for an example
  of a class derived from this one
  @short abstract interface to weighting schemes when codons differ at 2 positions
*/

/*!
  \fn Sequence::WeightingScheme2::Calculate(const std::string &codon1, const std::string &codon2)const =0;
  Calculate actually calculates the weights for each branch
  \param codon1 a std::string of length 3 representing a sense codon
  \param codon2 a std::string of length 3 representing a sense codon
*/

/*!
  \fn Sequence::WeightingScheme2::weights(void) const = 0;
  \return a double * of size 6 (1 value for each branch)
*/

/*!
    \class Sequence::WeightingScheme3 Sequence/WeightingSchemes.hpp
    \ingroup weights
    \warning Before deriving classes from this base, please study the implementations
    of Sequence::ThreeSubs.  You should also study Sequence::GranthamWeights3 for an example
    of a class derived from this one
    @short abstract interface to weighting schemes when codons differ at 3 positions
  */

/*!
  \fn Sequence::WeightingScheme3::Calculate(const std::string &codon1, const std::string &codon2)const =0;
  Calculate actually calculates the weights for each branch
  \param codon1 a std::string of length 3 representing a sense codon
  \param codon2 a std::string of length 3 representing a sense codon
*/

/*!
  \fn Sequence::WeightingScheme3::weights(void) const = 0;
  \return a double * of size 6 (1 value for each branch)
*/
namespace Sequence
  {


  class WeightingScheme2
    {
    private:
    public:
      explicit WeightingScheme2(void)
      {}
      virtual ~WeightingScheme2(void)
      {}
      virtual void Calculate(const std::string &codon1, const std::string &codon2) const =0;
      virtual double *weights(void) const = 0;
    };

  class WeightingScheme3
    {
    private:
    public:
      explicit WeightingScheme3(void)
      {}
      virtual ~WeightingScheme3(void)
      {}
      virtual void Calculate(const std::string &codon1, const std::string &codon2)const =0;
      virtual double *weights(void) const = 0;
    };
}
#endif
