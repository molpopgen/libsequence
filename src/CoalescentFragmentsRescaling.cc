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

#include <Sequence/Coalescent/FragmentsRescaling.hpp>
#include <Sequence/SimData.hpp>
#include <cassert>
#include <numeric>

namespace Sequence
{
  int sample_length( const std::vector< std::pair<int,int> > & fragments )
  /*!
    \return The sum of fragments[i].second for i=0 to i=fragments.size()-1
    \brief When simulating partially linked regions, return the total length of sample material that we are simulating
    \ingroup coalescent
  */
  {
    int sum = 0;
    for( std::vector< std::pair<int,int> >::const_iterator i = fragments.begin() ;
	 i < fragments.end() ; ++i )
      {
	sum += i->second;
      }
    return sum;
  }

  int total_length( const std::vector< std::pair<int,int> > & fragments )
  /*!
    \return The sum of fragments[i].first + fragments[i].second for i=0 to i=fragments.size()-1
    \brief When simulating partially linked regions, return the total length of the region
    \ingroup coalescent
  */
  {
    int sum = 0;
    for( std::vector< std::pair<int,int> >::const_iterator i = fragments.begin() ;
	 i < fragments.end() ; ++i )
      {
	sum += (i->first+i->second);
      }
    return sum;
  }

  void calculate_scales(const std::vector< std::pair<int,int> > & fragments,
			std::vector< std::pair<double,double> > * sample_scale,
			std::vector< std::pair<double,double> > * mutation_scale )
  /*!
    \brief This is a helper function that rescales physical distance in base pairs to
    continuous distance on the interval 0,1.  
    \param fragments  A vector of pairs, representing physical distance in bp.
    For each pair, the first element is the distance to the next fragment, and the
    second element is the length of the fragment.  For example, two 1kb fragments
    separated by 10kb would be represented by the pairs (0,1000) (10000,1000).
    \param sample_scale  This vector will be filled with values representing the
    positions of the fragments on the continuous interval, without any space betwen them.
    This is because we will actually do the simulation using a non-uniform genetic map
    to represent the high recombination rates between fragments
    \param mutation_scale This is a direct mapping of the data contained in \a fragments
    to the continuous scale, and can be used to rescale the positions of mutations
    \ingroup coalescent
  */
  {
    sample_scale->resize(fragments.size());
    mutation_scale->resize(fragments.size());
    //1. need to get the total length of the region simulated
    int ttl_len = total_length(fragments);
    int ttl_sample_len = sample_length(fragments);

    //2. need to convert the positions in fragments into scaled units
    //we'll store this as a pair of [beg,end) rather than distance, length
    double dummy = 0.;
    int ttl = 0;
    std::vector< std::pair<double,double> >::iterator si = sample_scale->begin(),
      mi = mutation_scale->begin();
    for( std::vector<std::pair<int,int> >::const_iterator i = fragments.begin() ;
	 i < fragments.end() ; ++i,++si,++mi )
      {
	si->first = dummy;
	si->second = double(i->second)/double(ttl_sample_len);
	dummy += double(i->second)/double(ttl_sample_len);
	mi->first = double(ttl)/double(ttl_len) + double(i->first)/double(ttl_len);
	mi->second = double(i->second)/double(ttl_len);
	ttl += i->first + i->second;
      }
  }

  void rescale_mutation_positions(SimData * d,
				  const std::vector< std::pair<double,double> > & sample_scale, 
				  const std::vector< std::pair<double,double> > & mutation_scale )
  /*!
    \brief Rescales the positions of the mutations in \a d from the scale given in \a sample_scale to that
    given in \a mutation_scale
    \note See documentation for calcualate_scales
    \ingroup coalescent
  */
  {
    assert(mutation_scale.size()==sample_scale.size());
    typedef std::vector<std::pair<double,double> >::const_iterator ci;
    for(SimData::pos_iterator p = d->pbegin() ; p < d->pend() ; ++p)
      {
	for( ci ss=sample_scale.begin(),ms=mutation_scale.begin() ; 
	     ss < sample_scale.end() && ms < mutation_scale.end() ; ++ss,++ms )
	  {
	    if( *p >= ss->first &&
		*p < (ss->first+ss->second) )
	      {
		//is in fragment i
		double delta = (*p-ss->first)/(ss->second);
		*p = ms->first + delta*ms->second;
		break;
	      }
	  }
      }
  }

  void rescale_arg( arg * sample_history,
		    const std::vector< std::pair<int,int> > & fragments )
  /*!
    \brief Rescales the beginnings of marginal trees in an ancestral recombination graph 
    from a genetic scale to a physical scale
    \ingroup coalescent
  */
  {
    for(arg::iterator ai = sample_history->begin() ; ai != sample_history->end() ; ++ai)
      {
	int ttl_len = 0;
	int ttl_sample_len=0;
	for( std::vector< std::pair<int,int> >::const_iterator f = fragments.begin() ;
	     f < fragments.end() ; ++f )
	  {
	    ttl_len += f->first;
	    if( ai->beg >= ttl_sample_len && ai->beg < (ttl_sample_len+f->second) ) 
	      //the i-th marginal is in fragment f
	      {
		double delta = double(ai->beg-ttl_sample_len)/double(f->second);
		ai->beg = ttl_len + int(delta*double(f->second));
		break;
	      }
	    else
	      {
		ttl_sample_len += f->second;
		ttl_len+=f->second;
	      }
	  }
      }
  }

  double integrate_genetic_map( const std::vector<chromosome> & sample,
				const int & current_nsam,
				const std::vector<double> & genetic_map,
				std::vector<double> * reclens)
  /*!
    \brief When simulating non-uniform recombination rates, the probability of recombination
    at each point in the simulation needs to be obtained by integrating over the genetic map and 
    the current sample configuration.  This function does that.
    \param sample the vector containing the current state of all chromosomes in the sample
    \param current_nsam the current sample size in the simulation
    \param genetic_map a vector containing rho/"link" for each link in the sample.  
    For the i-th base-pair in the chromosome, the "link" is the "space between" positions i and i+1.  
    The value of genetic_map[i] is therefore 4Nr between site i and i+1 (sometimes called 4Nr/"site").
    \param reclens a vector of doubles.  This vector will be resized to \a current_nsam in this function,
    and filled with \a current_nsam values, each of which is the sum(genetic_map[beg],genetic_map[end-1]) 
    for each chromosome in the sample, where beg and end are the first and last positions in each chromosome.  
    These data are needed by the function pick_spot (Sequence/Coalescent/Recombination.hpp).
    \return the cummulative recombination rate in the sample, which is obtained by integrating over 
    the ancestral material in the sample and the genetic map.
    \ingroup coalescent
  */
  {
    assert(current_nsam > 0);
    assert(!genetic_map.empty());
    reclens->resize(current_nsam);
    std::vector<double>::iterator ri = reclens->begin();
    double rrec=0.;
    for(std::vector<chromosome>::const_iterator chrom = sample.begin() ;
	chrom < (sample.begin()+current_nsam) ; ++chrom,++ri)
      {
	int beg = chrom->first(), end = chrom->last();
	assert(beg<=end);
	double rho_chrom = std::accumulate(genetic_map.begin()+beg,genetic_map.begin()+end,0.);
	*ri = rho_chrom;
	rrec += rho_chrom;
      }
    return rrec;
  }
} //ns Sequence
