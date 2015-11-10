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

#ifndef NDEBUG
#include <cmath>
#include <cfloat>
#endif
#include <cassert>
#include <algorithm>
#include <cctype>
#include <Sequence/Seq.hpp>
#include <Sequence/Translate.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/Comparisons.hpp>
#include <Sequence/RedundancyCom95.hpp>
#include <Sequence/Sites.hpp>
//divergence statistics for a pair of sequences

namespace Sequence
{
  struct SitesImpl
  {
    double _L0;
    double _L2S;
    double _L2V;
    double _L4;
    size_t seqlen;
    int maxhits;
    void siteinc (const RedundancyCom95 * sitesObj,
		  const std::string & codon1,const  std::string &codon2);
    void count_sites (const Sequence::Seq * sequence1,
		      const Sequence::Seq * sequence2,
		      const RedundancyCom95 * sitesObj);
    SitesImpl(const Sequence::RedundancyCom95 * sitesObj, const Sequence::Seq * seq1,
	      const Sequence::Seq * seq2, const int max) :
      _L0(0.),_L2S(0.),_L2V(0.0),_L4(0.0),
      seqlen(seq1->size()),
      maxhits(max)
    {
    }
  };
  
  Sites::Sites (const Sequence::RedundancyCom95 * sitesObj, const Sequence::Seq * seq1,
		const Sequence::Seq * seq2, const int max) :
    impl(std::unique_ptr<SitesImpl>(new SitesImpl(sitesObj,seq1,seq2,max)))
    /*!
      \param sitesObj an initialized object of type Sequence::RedundancyCom95
      \param seq1 a Sequence::Seq
      \param seq2 a Sequence::Seq
      \param max max number of substitutions per codon to analyze
      \note sequences must be of same length, this is checked by assert()
      \note sequence lengths must be multiples of 3, this is checked by assert()
    */
  {
  }

  Sites::~Sites (void)
  {}

  void
  SitesImpl::count_sites (const Sequence::Seq * sequence1,
			  const Sequence::Seq * sequence2,
			  const RedundancyCom95 * sitesObj)
  {				//now, it is correct
    size_t i = 0, j = 0;
    std::string codon1,codon2;
    codon1.resize(3);
    codon2.resize(3);
    for (i = 0; i <= seqlen - 3; i += 3)
      {
	for (j = 0; j <= 2; ++j)
	  {
	    codon1[j] = char(std::toupper((*sequence1)[i + j]));
	    codon2[j] = char(std::toupper((*sequence2)[i + j]));
	  }
	//We won't check the r.v. of NumDiffs here b/c we have made the 2 seqs the
	//same lenght.
	decltype(NumDiffs(codon1,codon2)) nc = NumDiffs (codon1, codon2);

	if (nc == 0)	//still need to count if there are 0 changes
	  siteinc (sitesObj,codon1, codon2);
	else if (maxhits <= 3 && nc == 1)
	  siteinc (sitesObj,codon1, codon2);
	else if (maxhits == 2 && nc <= 2)
	  siteinc (sitesObj,codon1, codon2);
	else if (maxhits == 3 && nc <= 3)
	  siteinc (sitesObj,codon1, codon2);
      }
  }

  void
  SitesImpl::siteinc (const RedundancyCom95 * sitesObj,
		      const std::string & codon1,
		      const std::string & codon2)
  {
    //skip ambiguous

    if ( std::find_if(codon1.begin(),codon1.end(),ambiguousNucleotide()) == codon1.end() &&
	 std::find_if(codon2.begin(),codon2.end(),ambiguousNucleotide()) == codon2.end() )
      {
	_L0 += (sitesObj->L0_vals(codon1) + sitesObj->L0_vals(codon2))/2.0;
	_L2S += (sitesObj->L2S_vals(codon1) + sitesObj->L2S_vals(codon2))/2.0;
	_L2V += (sitesObj->L2V_vals(codon1) + sitesObj->L2V_vals(codon2))/2.0;
	_L4 += (sitesObj->L4_vals(codon1) + sitesObj->L4_vals(codon2))/2.0;
      }
  }
	 

  double Sites::L0(void) const
  /*!
    \return alignment length in terms of non-degenerate sites
  */
  {
    return impl->_L0;
  }
  double Sites::L2S(void) const
  /*!
    \return alignment length in terms of transitional-degenerate sites
  */
  {
    return  impl->_L2S;
  }
  double Sites::L2V(void) const
  /*!
    \return alignment length in terms of transversional-degenerate sites
  */
  {
    return impl->_L2V;
  }
  double Sites::L4(void) const
  /*!
    \return alignment length in terms of fourfold-degenerate sites
  */
  {
    return impl->_L4;
  }
	 
}
