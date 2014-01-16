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

//Kimura80.cc Copyright 2001 Kevin Thornton, k-thornton@uchicago.edu
#include <cmath>
#include <cfloat>
#include <cctype>
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#include <Sequence/Seq.hpp>
#include <Sequence/Comparisons.hpp>
#include <Sequence/SeqEnums.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/Kimura80.hpp>

/*!
  \defgroup divergence Divergence Statistics
*/

//using std::isnan;
using std::log;
using std::isfinite;

namespace Sequence
  {
  Kimura80::Kimura80 (const Sequence::Seq * seqa, const Sequence::Seq * seqb):
    seqlen (seqa->length())
      /*!
        \param seqa an object of type Sequence::Seq
        \param seqb an object of type Sequence::Seq
        \exception Sequence::SeqException if sequences are of different lengths
      */
  {
    if (seqa->length () != seqb->length ())
      throw SeqException ("Sequence::Kimura80::Kimura80(): constructor called with two sequence objects of unequal length");
    num_Ts = 0;
    num_Tv = 0;
    divergence = 0.0;
    P = 0.0;
    Q = 0.0;
    sites_compared = 0;
    Compute (seqa,seqb);
  }

  void Kimura80::Compute (const Sequence::Seq *seq1, const Sequence::Seq *seq2)
  {
    unsigned i;

    unsigned ungapped_sites = 0;
    for (i = 0; i < seqlen; ++i)	//iterate over the sequence
      {
        int type = 0;
        if (NotAGap((*seq1)[i]) && NotAGap((*seq2)[i]))
          {
            ++ungapped_sites;
            if (std::toupper((*seq1)[i]) != 
		std::toupper((*seq2)[i]))	//if the sites differ at that position
              {
                type = TsTv ((*seq1)[i], (*seq2)[i]);	//check if difference is Ts or Tv
              
		if (type == Mutations(Ts))
		  {	//Ts
		    ++num_Ts;
		  }
		else if (type == Mutations(Tv))
		  {	//Tv
		    ++num_Tv;
		  }
	      }
          }
      }
    //make sure we use the right denominator

    sites_compared = (ungapped_sites < seqlen) ? ungapped_sites : seqlen;
    //P and Q are the proportions of Ts and Tv changes observed
    P = double (num_Ts) / double (sites_compared);
    Q = double (num_Tv) / double (sites_compared);

    //Kimura's formula
    if (fabs(1.0 - 2.0 * P - Q) > DBL_EPSILON)
      {
        divergence = -1.0 * 0.5 * log ((1.0 - 2.0 * P - Q)
                                       *  pow ((1 - 2.0 * Q), 0.5));
      }
    else
      {
        divergence = 0.;
      }
    //a correction for extremely low observed values
    if (divergence <= 0.0-DBL_EPSILON)
      divergence = 0.0;
  }

  double
  Kimura80::K (void)
  /*!
    \return the distance between the two sequences.
    \note 999.0 is returned as a warning value.
    This can be necessary if sites are saturated,
    which implies that divergence cannot be calculated
  */
  {
    if (!isfinite(divergence))
      return (999.0);
    return (divergence);
  }

    size_t
    Kimura80::sites (void)
    /*!
      \return the number of sites compared, excluding gaps, missing data, etc.
    */
    {
      return (sites_compared);
    }
}
