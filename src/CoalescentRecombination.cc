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

#include <Sequence/Coalescent/Recombination.hpp>
#include <cassert>

#ifndef NDEBUG
namespace 
{
  bool arg_is_sorted( const Sequence::arg * sample_history )
  {
    for( Sequence::arg::const_iterator i = sample_history->begin() ;
	 i != sample_history->end() ;
	 ++i )
      {
	Sequence::arg::const_iterator j=i;
	j++;
	if( j != sample_history->end() )
	  {
	    if(i->beg > j->beg)
	      return false;
	  }
      }
    return true;
  }
}
#endif

namespace Sequence
{
  std::pair<int,int> pick_uniform_spot(const double & random_01,
				       const int & nlinks,
				       std::vector<chromosome>::const_iterator sample_begin,
				       const unsigned & current_nsam)
  /*!
    @brief Pick a crossover point for the model where recombination rates are constant
    across a recion.
    Picks a positions uniformly amongst all chromosomes at which a recombination event
    will occur.
    \param random_01 a random uniform deviate U[0,1)
    \param nlinks the number of links currently in the simulation
    \param sample_begin an iterator pointing to the beginning of the sample
    \param current_nsam the current sample size in the simulation
    \return a pair of integers containing the index of the recombinant chromosome 
    (.first), and the position at which the crossover will occur (.second)
    \ingroup coalescent
  */
  {
    int pos = int(random_01*double(nlinks))+1;
    int recombinant = 0,len=0;
    while( (sample_begin+recombinant) <
	   (sample_begin+current_nsam) )
      {
	len = (sample_begin+recombinant)->links();
	if(pos <= len)break;
	pos -= len;
	recombinant++;
      }
    int rpos = (sample_begin+recombinant)->begin()->beg + pos - 1;
    return std::make_pair(recombinant,rpos);
  }

  int crossover( const int & current_nsam,
		 const int & chromo,
		 const int & pos,
		 std::vector<chromosome> * sample,
		 arg * sample_history)
  /*!
    @brief Recombination function.
    \param current_nsam the current sample size in the simulation
    \param chromo the chromosome on which the crossover event is to occur
    \param pos the crossover event happens between sites pos and pos+1 (0<= pos < nsites)
    \param sample the sample of chromosomes being simulated
    \param sample_history the genealogy of the sample
    \return the number of links lost due to the crossover event
    \note as the type arg is based on std::list, and insertions into lists are done
    in constant time, this routine keeps the ancestral recombination graph sorted
    \ingroup coalescent
  */
  {
    std::vector<chromosome>::iterator sbegin = sample->begin();

    //1. Is pos within a segment, or between segments?
    bool within = false;
    chromosome::iterator seg = (sbegin+chromo)->begin();
    for( ; pos >= seg->end ;++seg){};
    assert(seg != (sbegin+chromo)->end());
    within = (pos>=seg->beg) ? true:false;

    //2. make new chromosome segments for right-hand end
    size_t ns = (sbegin+chromo)->nsegs - (seg-(sbegin+chromo)->begin());
    std::vector<segment> rtsegs(seg,seg+ns);

    //3. edit vector of segments for left-hand end
    (sbegin+chromo)->nsegs = (int((seg)-(sbegin+chromo)->begin())) + int(within);
    assert( (sbegin+chromo)->nsegs > 0 );

    //4. make sure begs and ends are happy
    if(within)
      {
	rtsegs.begin()->beg = pos+1;
	seg->end = pos;
      }
    else
      {
	rtsegs.begin()->beg = seg->beg;
      }

    //rv is the number of links lost due to the crossover event
    int rv = rtsegs.begin()->beg - ((sbegin+chromo)->end()-1)->end;

    //5. insert a new chromosom into the sample
    //WARNING: all pointers and iterators declared above
    //should be considered invalidated!
    int tpop=(sbegin+chromo)->pop;
    sample->insert(sbegin+current_nsam,chromosome(rtsegs,tpop));
    //code below causes too much RAM usage  (not sure why...)
//     if(std::vector<chromosome>::size_type(current_nsam)+1 > sample->size())
//       {
//         sample->resize(current_nsam+100);
//         sbegin = sample->begin();
//       }
//     *(sbegin+current_nsam) = chromosome(rtsegs,tpop);
    assert( (sample->begin()+current_nsam)->nsegs == rtsegs.size() 
	    && rtsegs.size() == ns );

    //6. insert a new marginal tree if necessary
    if(within == true)
      {
	int beg_new_marg = rtsegs.begin()->beg;
	
	//find place in arg that is affected
	arg::iterator argbeg = sample_history->begin();
	arg::iterator titr=argbeg;
	titr++;
	for( ; titr != sample_history->end()
	       && beg_new_marg > titr->beg-1 ; argbeg++,titr++ ){};
	assert(argbeg!=sample_history->end());
	if(argbeg->beg != beg_new_marg)
	  {
	    arg::iterator argt = sample_history->insert(titr,*(argbeg));
	    argt->beg = beg_new_marg;
	    assert(arg_is_sorted(sample_history));
	  }
      }
    return rv;
  }
}
