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

#ifndef REDCOM95_H
#define REDCOM95_H
/*! \file RedundancyCom95.hpp
  @brief  class Sequence::RedundancyCom95
*/

/*!
  \class Sequence::RedundancyCom95 Sequence/RedundancyCom95.hpp
  \ingroup kaks
  This class performs "on-the-fly" tabulations of the redundancy of a genetic code,
  in terms of whether or not a transition or transversion at any site in the codon
  changes the amino acid encoded by that codon.  It supports
  whatever codes libsequence supports, because it uses Sequence::Translate to determine whether
  a mutation at a certain position in a codon will encode the same amino acid or not.\n
  \n
  This should be considered a "buried" class, in the sense that it is, properly speaking, an 
  implementation detail specific to Sequence::Comeron95. Another possible use for this class
  (which I will implement), is for analysis of silent/replacement polymorphism in
  population-genetic  data.
 
  \note any variable of type std::string that is named codon implicity assumes codon.length()==3
  @short Calculate redundancy of a genetic code using Comeron's counting scheme
*/
#include <string>
#include <memory>
#include <Sequence/SeqEnums.hpp>
#include <Sequence/SeqExceptions.hpp>

using Sequence::GeneticCodes;

namespace Sequence
{
  class  RedundancyCom95impl;
  class RedundancyCom95
  {
  private:
    std::auto_ptr<RedundancyCom95impl> impl;
  public:
    explicit RedundancyCom95 (const Sequence::GeneticCodes &genetic_code = Sequence::UNIVERSAL);
    ~RedundancyCom95(void);
    //counting routines, return values from private matrices of same name
    double FirstNon (const std::string &codon) const ;
    double First2S (const std::string &codon) const ;
    double First2V (const std::string &codon) const ;
    double ThirdNon (const std::string &codon) const ;
    double ThirdFour (const std::string &codon) const ;
    double Third2S (const std::string &codon) const ;
    double Third2V (const std::string &codon) const ;
    double L0_vals (const std::string &codon) const ;
    double L2S_vals (const std::string &codon) const ;
    double L2V_vals (const std::string &codon) const ;
    double L4_vals (const std::string &codon) const ;
  };
}

#endif
