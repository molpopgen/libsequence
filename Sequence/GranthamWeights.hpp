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

#ifndef __GRANTHAMWEIGHTS_H__
#define __GRANTHAMWEIGHTS_H__

/*! \file GranthamWeights.hpp
  @brief declaration of classes to weight codons by Grantham distance (i.e. for Sequence::Comeron95). Declares 
  Sequence::GranthamWeights2 and  Sequence::GranthamWeights3
*/

/*!
  \class Sequence::GranthamWeights2 Sequence/GranthamWeights.hpp
  \ingroup weights
  @short Weights paths by Grantham's distances for codons differing at 2 sites
*/

/*!
  \class Sequence::GranthamWeights3 Sequence/GranthamWeights.hpp
  \ingroup weights
  @short Weights paths by Grantham's distances for codons differing at 3 sites
*/
#include <Sequence/SeqEnums.hpp>
#include <Sequence/WeightingSchemes.hpp>

namespace Sequence
  {
  class Grantham;
  class GranthamWeights2 : public WeightingScheme2
    {
    private:
      Sequence::GeneticCodes code;
      mutable double __weights[2];//logical const
    public:
      explicit GranthamWeights2(Sequence::GeneticCodes genetic_code = Sequence::UNIVERSAL);
      ~GranthamWeights2(void);
      void Calculate(const std::string &codon1, const std::string &codon2) const;
      double *weights(void) const;
    };

  class GranthamWeights3 : public WeightingScheme3
    {
    private:
      GeneticCodes code;
      mutable double __weights[6];//logical const
    public:
      explicit GranthamWeights3(Sequence::GeneticCodes genetic_code = Sequence::UNIVERSAL);
      ~GranthamWeights3(void);
      void Calculate(const std::string &codon1, const std::string &codon2) const;
      double *weights(void) const;
    };
}
#endif
