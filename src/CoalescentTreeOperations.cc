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

#include <Sequence/Coalescent/TreeOperations.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <cassert>
#include <cmath>
#include <limits>

namespace Sequence
{
  double total_time( const marginal::const_iterator beg,
		     const int & nsam )
  /*!
    @brief Calculate total time on a marginal tree
    \param beg A pointer to the beginning of a marginal tree,
    i.e. the return value of marginal::begin()
    \param nsam the total sample size simulated
    \return The total time on the tree.
    \note The scaling of time in the simulation is up to you
    \ingroup coalescent
  */
  {
    double t = (beg+2*nsam-2)->time;
    int i = nsam;
    while( i < 2*nsam-1 )
      {
	t += (beg+i)->time;
	++i;
      }
    return t;
  }
  
  int pick_branch(marginal::const_iterator beg,
		  const int & nsam,
		  const double & rtime)
  /*!
    @brief pick a random branch of a marginal tree
    \param beg A pointer to the beginning of a marginal tree,
    i.e. the return value of marginal::begin()
    \param nsam the total sample size simulated
    \param rtime a (preferably random) double between 0 and the
    total_time on the marginal tree from which \a beg is the iterator
    \ingroup coalescent
  */
  {
    double time=0.;
    int i=0;
    //we iterate over all-but-the-last branch
    //this is a guard in case rtime == ttime,
    //which can cause numerical issues in 
    //the internal if statement
    for(i=0;i<2*nsam-3;i++)
      {
	time += (beg + (beg+i)->abv)->time - (beg+i)->time;
	if(time >= rtime) return i;
      }
    return i;
  }

  std::vector<int> get_all_descendants (marginal::const_iterator beg,
					const int & nsam,
					const int & branch)
  /*!
    @brief Find all the descendants of a branch on a marginal tree
    \param beg A pointer to the beginning of a marginal tree,
    i.e. the return value of marginal::begin()
    \param nsam the total sample size simulated
    \param branch the index of the branch of the tree whose descendants you want.
    \note branch must be <= 2*nsam-2, which is checked by assert
    \ingroup coalescent
  */
  {
    assert( branch <= 2*nsam-2 );
    std::vector<int> desc;
    int j;
    for(int i=0;i<nsam;++i)
      {
	for(j=i;j<branch;j=(beg+j)->abv){};
	if(j==branch) desc.push_back(i);
      }
    return desc;
  }

  bool is_descendant( marginal::const_iterator beg,
		      const int & ind,
		      const int & branch )
  /*!
    @brief Ask if a tip of a tree is a descendant of a particular branch
    \param beg A pointer to the beginning of a marginal tree,
    i.e. the return value of marginal::begin()
    \param ind the index of the putative descendant node
    \param branch the index of the branch of the tree which may be the ancestor
    of \a ind
    \note This function does not check whether ind or branch go out of bounds,
    and so the programmer must ensure that both values are <= 2*nsam-2, where
    nsam is the total sample size simulated
    \ingroup coalescent
  */
  {
    int i=ind;
    while(i<branch)
      {
	assert( (beg+i)->abv != -1 );
	i = (beg+i)->abv;
      }
    return i==branch;
  }

  double total_time_on_arg( const Sequence::arg & sample_history,
			    const int & total_number_of_sites )
  /*!
    @brief Returns the total time on an ancestral recombination graph.
    \param sample_history an ancestral recombination graph
    \param total_number_of_sites the number of "sites" simulated on the ARG
    \return The total time on an ancestral recombination graph.
    \note The time is in terms of whatever units are recorded on the nodes of the mariginals of the ARG
    \throws Sequence::SeqException if the beginning of any marginal tree is >= total_number_of_sites
    \ingroup coalescent
   */
  {
    
    size_t seg=0,nsegs=sample_history.size();
    Sequence::arg::const_iterator i = sample_history.begin(),j=i;
    j++;
    
    double tt = 0.;
    
    for( ; seg < nsegs ; ++seg,++i,++j )
      {
	//check that things are sane
	assert( i->beg < total_number_of_sites );
	if( i->beg > total_number_of_sites )
	  {
	    throw Sequence::SeqException("Sequence::total_time_on_arg - beginning of marginal tree > number of sites simulated");
	  }
	//get the beginning and end of the current marginal tree
	int end = (seg<nsegs-1) ? j->beg : total_number_of_sites;
	int beg = i->beg;
	
	//increment total time on tree as a weighted average
	//of the times on all the marginals
	tt += Sequence::total_time(i->begin(),i->nsam) * 
	  double(end-beg)/double(total_number_of_sites);
      }
    return tt;
  }

  void minimize_arg( Sequence::arg * sample_history )
  /*!
    Takes an arg (Ancestral Recombination Graph) and removes redundant marginal trees.
    Specifically, for two adjacent trees i and j (j->beg > i->beg), j is removed from the arg
    if the topology and branch lengths of i and j are identical.
    \param sample_history the arg to minimize
  */
  {
    arg::iterator i = sample_history->begin(),
      j=i;
    ++j;
    assert(j->beg > i->beg);
    size_t nsegs = sample_history->size();
    for(size_t seg=0;seg<nsegs-1;++seg,++i,++j)
      {
	bool same_tree = true;
	for( marginal::iterator mi = i->end()-1,
	       mj=j->end()-1 ;
	     same_tree == true && (mi >= i->begin() && mj >= j->begin()) ;
	     --mi,--mj )
	  {
	    if(std::fabs( mi->time - mj->time ) > std::numeric_limits<double>::epsilon() ||
	       mi->abv != mj->abv )
	      {
		same_tree = false;
	      }
	  }
	if(same_tree == true)
	  {
	    sample_history->erase(j);
	    nsegs--;
	    j=i;
	    i--;
	    seg--;
	}
      }
  }

#ifndef DOXYGEN_SKIP
  struct sfs_times_impl
  {
    std::vector<double> times;
    double tt;
    size_t nbins;
    sfs_times_impl() : times(std::vector<double>()),tt(0.),nbins(0){}
    sfs_times_impl( const sfs_times_impl & s ) :
      times(s.times),tt(s.tt),nbins(s.nbins) {}
    sfs_times_impl(arg::const_iterator & sample_history_beg,
		   const arg::size_type & nsegs,
		   const int & total_nsites_simulated,
		   bool folded);
    inline bool operator==(const sfs_times_impl & rhs) const
    {
      return (this->times==rhs.times &&
	      std::fabs(this->tt-rhs.tt)<=std::numeric_limits<double>::epsilon() &&
	      this->nbins == rhs.nbins);
    }
  };

  sfs_times_impl::sfs_times_impl(arg::const_iterator & sample_history_beg,
				 const arg::size_type & nsegs,
				 const int & total_nsites_simulated,
				 bool folded)
    : times(std::vector<double>(((folded==false) ? sample_history_beg->nsam-1 : 
				 (sample_history_beg->nsam/2)),0.)),
	    tt(0.),
	    nbins(((folded==false) ? sample_history_beg->nsam-1 : 
		   (sample_history_beg->nsam/2)))
  {
    arg::const_iterator i = sample_history_beg,j=i;
    j++;
    int beg,end;
    double t;
    for(unsigned seg = 0; seg < nsegs ; ++seg,++i,++j)
      {
	end = (seg<nsegs-1) ? j->beg : total_nsites_simulated;
	beg = i->beg;
	const double scale = double(end-beg)/double(total_nsites_simulated);
	//iterate over tips
	const marginal::const_iterator treebeg = i->begin();
	for(int tip = 0 ; tip < i->nsam ; ++tip)
	  {
	    //tips lead to singletons
	    t = ( (treebeg+((treebeg+tip )->abv))->time - 
		  (treebeg+tip)->time)*scale;
	    tt += t;
	    times[0] += t;
	  }
	//iterate over internal nodes;
	for(int node = i->nsam ; node < 2*(i->nsam)-2 ; ++node)
	  {
	    std::vector<int> descendants = get_all_descendants(treebeg,i->nsam,node);
	    t = ( (treebeg+((treebeg+node )->abv))->time - 
		  (treebeg+node)->time)*scale;
	    std::vector<double>::size_type index = (folded==false) ? descendants.size()-1 : 
	      std::min(descendants.size(),i->nsam-descendants.size())-1;
	    times[index] += t;
	    tt+=t;
	  }
      }
  }
#endif

  sfs_times::sfs_times() : 
    impl(std::auto_ptr<sfs_times_impl>(new sfs_times_impl()))
		      /*! Construct empty object */
  {
  }

  sfs_times::sfs_times(const sfs_times & sfst) : impl(new sfs_times_impl(*(sfst.impl)))
		      /*!
			Copy constructor
		      */
  {
  }

  sfs_times::sfs_times(arg::const_iterator sample_history_beg,
		       const arg::size_type & nsegs,
		       const int & total_nsites_simulated,
		       bool folded) : 
    impl(new sfs_times_impl(sample_history_beg,nsegs,total_nsites_simulated,folded))
		      /*!
			\param beg pointer to the first element in an Ancestral Recombination Graph (ARG)
			\param nsegs number of marginal trees in the ARG.  This will be equal to the return
			value of the member function arg::size()
			\param total_sites_simulated  This will equal 1 for models without recombination,
			else the nsites value used.
		       */
			
  {
  }

  sfs_times::~sfs_times()
  {
  }

  double sfs_times::operator[](std::vector<double>::size_type const & i) const
  /*!
    \return the total time leading to the i-th bin in the site frequency spectum
    \note the indexing for this object is from 1 <= x <= nsam-1, where nsam is the total
    sample size in the arg
  */
  {
    assert(i>0 && i<=impl->nbins);
    return impl->times[i-1];
  }

  sfs_times & sfs_times::operator=(const sfs_times & rhs)
  /*!
    assignment operator
  */
  {
    if(*this==rhs) return *this;
    this->impl->times = rhs.impl->times;
    this->impl->tt = rhs.impl->tt;
    this->impl->nbins = rhs.impl->nbins;
    return *this;
  }

  bool sfs_times::operator==(const sfs_times & rhs) const
  /*!
    \return *(this->impl) == *(rhs.impl)
  */
  {
    return ( *(this->impl) == *(rhs.impl) );
  }

  double sfs_times::ttime() const
  /*!
    \return total time in the ARG.  useful for checking...
  */
  {
    return impl->tt;
  }

  size_t sfs_times::size() const
  /*!
    \return the number of bins in the site frequency spectrum.  For a sample of size n,
    this is n-1 for the unfolded site frequency spectrum, and n/2 for the folded.
  */
  {
    return impl->nbins;
  }

  sfs_times::const_iterator sfs_times::begin() const
  /*!
    \return a const iterator to the first time (singletons)
  */
  {
    return impl->times.begin();
  }

  sfs_times::const_iterator sfs_times::end() const
  /*!
    \return a const iterator to 1 past the last time
  */
  {
    return impl->times.end();
  }
}
