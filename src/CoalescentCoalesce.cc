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

#include <Sequence/Coalescent/Coalesce.hpp>
namespace Sequence
{
  bool isseg( chromosome::const_iterator seg, const int nsegs,
	      const int pos, int * offset )
  /*!
    @brief ask if a chromosome beginning at seg and containing nsegs
    contains a segment containing the position pos
    \param seg a pointer to a segment of a chromosome (this should 
    be the 1st segment, such as the return value of chromosome::begin())
    \param nsegs the number of segs in the chromosome pointed to by \a seg
    \param offset a pointer to an integer.  This integer is used 
    for repeated pointer arithmetic, and should be initalized to 0
    before the first call.
    \param pos a position a long a chromosome.  This function
    asks if pos is contained in the ancestral material of the
    chromosome whose segments begin at \a seg
    \return true if a segment exists that contains the point pos
    \note only used by the function coalesce
    \ingroup coalescent
  */
  {
    for( ; (seg+*offset)->beg <= pos && (*offset) < nsegs; ++(*offset) )
      if ((seg+*offset)->end>=pos) return true;
    return false;
  }

  int coalesce(const double & time,
	       const int & ttl_nsam,
	       const int & current_nsam,
	       const int & c1,
	       const int & c2,
	       const int & nsites,
	       int * nlinks,
	       std::vector<chromosome> * sample,
	       arg * sample_history)
  /*!
    @brief Common ancestor routine for coalescent simulation.  Merges
    chromosome segments and updates marginal trees.

    Common ancestor routine for coalescent simulation. This routine performs
    the merging of two lineages by a coalescent event.  Such merges usually 
    require two sorts of operations.  The first is an update to the segments
    contained in a chromosome, and the second is an update of the nodes
    on a marginal tree.
    \param time the time at which the coalecent event is occuring
    \param ttl_nsam the total sample size being simulated
    \param current_nsam the current sample size in the simulation
    \param c1 the array index of the first chromosome involved in the coalescent event
    \param c2 the array index of the second chromosome involved in the coalescent event
    \param nsites the total mutational length of the region begin simulated.  In the
    language of Hudson (1983), this is the number of infinitely-many-alleles loci in the
    simulation.
    \param nlinks a pointer to the number of "links" currently in the simulation.
    A link is the region between two sites, such that a chromosome currently with k sites has 
    k-1 links
    \param sample a pointer to the vector of chromosomes which makes up the sample
    \param sample_history a pointer to the ancestral recombination graph
    \return the decrease in current_nsam due to the coalescent event.  Usually, 
    the return value is 1.  Sometimes, however, it is two, when the two chromosomes
    being merged have no ancestral material on the same marginal tree.
    \ingroup coalescent
  */
  {
    bool yes1,yes2;
    int ch1=(c1<c2)?c1:c2, ch2=(c2>c1)?c2:c1;

    assert( (sample->begin()+ch1)->nsegs>0 );
    assert( (sample->begin()+ch2)->nsegs>0 );

    std::vector<chromosome>::iterator sbegin=sample->begin();

    chromosome::iterator ch1beg = (sbegin+ch1)->begin(),
      ch2beg=(sbegin+ch2)->begin();
    int seg1=0,seg2=0;

    segment * tsp = (segment *)malloc(sample_history->size()*sizeof(segment));
    int tseg = -1;

    //iterate over marginal histories
    int k=0,nsegs=int(sample_history->size());
    arg::iterator imarg = sample_history->begin(),
      jmarg=imarg;
    jmarg++;
    for ( ; k<nsegs ; ++imarg,++jmarg,++k )
      {
	//ask if chromosomes ch1 and ch2 have segments
	//that are part of the i-th marginal history
	yes1 = isseg(ch1beg,(sbegin+ch1)->nsegs,imarg->beg,&seg1);
	yes2 = isseg(ch2beg,(sbegin+ch2)->nsegs,imarg->beg,&seg2);
	if( yes1 || yes2 )
	  {
	    tseg++;
	    (tsp+tseg)->beg = imarg->beg;
	    (tsp+tseg)->end = (k<(nsegs-1)) ? jmarg->beg-1 : nsites-1;

	    if(yes1 && yes2)
	      {
		imarg->nnodes++;
		if( imarg->nnodes >= (2*ttl_nsam-2))
		  {
		    tseg--;
		  }
		else
		  {
		    (tsp+tseg)->desc = imarg->nnodes;
		  }
		marginal::iterator mi = imarg->begin();
		(mi+imarg->nnodes)->time = time;
		(mi+(ch1beg+seg1)->desc)->abv = imarg->nnodes;
		(mi+(ch2beg+seg2)->desc)->abv = imarg->nnodes;
		assert( (mi+(ch1beg+seg1)->desc)->abv <= int(2*ttl_nsam-2) );
		assert( (mi+(ch2beg+seg2)->desc)->abv <= int(2*ttl_nsam-2) );
	      }
	    else
	      {
		(tsp+tseg)->desc = (yes1==true) ? (ch1beg+seg1)->desc : (ch2beg+seg2)->desc;
		assert( (tsp+tseg)->desc < 2*ttl_nsam-1 );
	      }
	  }
      }
    *nlinks -= (sbegin+ch1)->links();
    int flag=0;
    if(tseg < 0)
      {
	free(tsp);
	(sbegin+ch1)->swap_with(*(sbegin+current_nsam-1));
	if(ch2 == current_nsam-1)
	  {
	    ch2=ch1;
	  }
	flag=1;
	assert( (sbegin+ch1)->nsegs>0 );
	assert( (sbegin+ch2)->nsegs>0 );
      }
    else
      {
	assert( (sbegin+ch1) < sample->end() );
	(sbegin+ch1)->assign_allocated_segs(tsp,tseg+1);
	*nlinks += (sbegin+ch1)->links();
      }

    *nlinks -= (sbegin+ch2)->links();
    (sbegin+ch2)->swap_with(*(sbegin+current_nsam-1-flag));
    return ((tseg<0)?2:1);  
  }
}
