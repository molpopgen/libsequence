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

/*! \file Kimura80.hpp
  @brief declaration of Sequence::Kimura80
*/

/*! \class Sequence::Kimura80 Sequence/Kimura80.hpp
  \ingroup divergence
  Calculate a measure of sequence divergence using Kimura's 1980 method.\n
  The reference is: Kimura, M (1980) J. Mol. Evol 16: 111-120.\n
  The calculation only depends on 3 numbers:\n
  1.) the number of sites in the sequence\n
  2.) the number of transitions between the two sequences\n
  3.) the number of transversions between the two sequences\n
  \n
  The implementation of this class does the following:\n
  1.) compare each position in both sequences, counting transitions and transversions\n
  2.) calculate distance using Kimura's formula\n
  \n
  \exception Sequence::SeqException if the two sequences are of unequal length.
 
  @short Kimura's 2-parameter distance
*/
#ifndef KIMURA80_H
#define KIMURA80_H

namespace Sequence
  {
  class Seq;
  class Kimura80
    {
    private:
      unsigned num_Ts, num_Tv;
      size_t seqlen;		//total sequence length
      size_t sites_compared;	//number of ungapped sites in the data
      void Compute (const Sequence::Seq *seq1, const Sequence::Seq *seq2);
      double divergence, P, Q;
    public:
      explicit Kimura80 (const Sequence::Seq * seqa,const  Sequence::Seq * seqb);
      double K ();
      size_t sites (void);
    };
}
#endif
