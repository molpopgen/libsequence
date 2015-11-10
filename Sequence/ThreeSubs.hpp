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

#ifndef THREESUBS_H
#define THREESUBS_H
/*! \file ThreeSubs.hpp
  @brief used by Sequence::Comeron95, class Sequence::ThreeSubs calculates divergence between codons that differ at three sites
*/
/*!
  \class Sequence::ThreeSubs Sequence/ThreeSubs.hpp
  \ingroup kaks
  \image html 3subs.jpg "Pathways when codons differ at 3 sites"
  \image latex 3subs.pdf "Pathways when codons differ at 3 sites" width=10cm
  A function object to obtain divergence statistics for Comeron's method for codons that differ at three positions. 
  Alternate paths are weighted by Grantham's distances. Used by Sequence::Comeron95
  \note You really need to study the above figure if you want to implement your own scheme to weight these
  pathways (done by deriving a class from Sequence::WeightingScheme3),
  you must know how this class names and assigns the intermediate steps and the branches.  In other 
  words, you must study the implementation of this class.  You may also want to refer to the
  implementations of Sequence::GranthamWeights3 for an example
  of implementing weighting schemes.
  @short Deal with codons differing at all 3 positions
*/
#include <Sequence/PathwayHelper.hpp>
#include <Sequence/WeightingSchemes.hpp>
#include <Sequence/SeqEnums.hpp>
#include <string>
#include <memory>
namespace Sequence
{
  class RedundancyCom95;
  struct ThreeSubsImpl;
  class ThreeSubs
  {
  private:
    std::unique_ptr<ThreeSubsImpl> impl;
  public:
    explicit ThreeSubs(void);
    void operator() (const RedundancyCom95 * sitesObj,
		     const std::string &codon1, const std::string &codon2,
		     const Sequence::WeightingScheme3 *weights3);
    ~ThreeSubs(void);
    //! \return number of transitions at non-degenerate sites in the codon
    double P0 (void) const;
    //! \return number of transitions at transitional-degenerate sites in the codon
    double P2S (void) const;
    //! \return number of transitions at transversional-degenerate sites in the codon
    double P2V (void) const;
    //! \return number of transitions at fourfold-degenerate sites in the codon
    double P4 (void) const;
    //! \return number of transversions at non-degenerate sites in the codon
    double Q0 (void) const;
    //! \return number of transversions at transitional-degenerate sites in the codon
    double Q2S (void) const;
    //! \return number of transversions at transversional-degenerate sites in the codon
    double Q2V (void) const;
    //! \return number of transversions at fourfold-degenerate sites in the codon
    double Q4 (void) const;
  };
}
#endif
