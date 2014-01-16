// Code for the -*- C++ -*- 
#ifndef __SEQUENCE_COALESCENT_BITS_MUTATION_TCC__
#define __SEQUENCE_COALESCENT_BITS_MUTATION_TCC__

#include <Sequence/Coalescent/TreeOperations.hpp>
#include <boost/ref.hpp>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>
#include <string>
#include <iostream>

namespace Sequence
{
#ifndef DOXYGEN_SKIP
  template<typename uniform_generator>
  void add_S_inf_sites_details( uniform_generator & uni ,
				marginal::const_iterator history,
				const double & tt,
				const int & beg, const int & end,
				const int & nsam, const int & nsites,
				const int & S ,
				const int & first_snp_index,
				gamete_storage_type * gametes )
  {
    assert(S>=0);
    assert(tt >= 0.);
    assert(nsam > 0);
    assert(nsites > 0);
    for(int snp = first_snp_index ; snp < first_snp_index+S ; ++snp)
      {
	double pos = uni(beg,end)/double(nsites);
	gametes->first[snp]=pos;
	int branch = pick_branch(history,nsam,uni(boost::cref(0.),tt));
	for(int ind=0;ind<nsam;++ind)
	  {
	    if( is_descendant(history,ind,branch) )
	      {
		gametes->second[ind][snp] = '1';
	      }
	    else
	      {
		gametes->second[ind][snp] = '0';
	      }
	  }
      }
  }

  template<typename poisson_generator,
	   typename uniform_generator>
  int infinite_sites_details( poisson_generator & poiss,
			      uniform_generator & uni,
			      gamete_storage_type * gametes,
			      const int & nsites,
			      const arg & history,
			      const double & theta)
  {
    assert(theta >= 0.);
    int snp_index_i=0;
    int ttlS=0;
    unsigned seg=0,nsegs=history.size();
    arg::const_iterator i = history.begin(),j=i;
    j++;
    for( ; seg < nsegs ; ++i,++j,++seg)
      {
	double tt = total_time(i->begin(),i->nsam);
	int end = (seg<nsegs-1) ? j->beg : nsites;
	int beg = i->beg;
	double mean = double(end-beg)*theta*tt/double(nsites);
	int S = poiss(mean);
	if(S>0)
	  {
	    ttlS+=S;
	    bool add=false;
	    while(ttlS >= MAX_SEGSITES)
	      {
		MAX_SEGSITES += MAX_SEGS_INC;
		add = true;
	      }
	    if(add)
	      {
		gametes->first.resize(MAX_SEGSITES);
		for(typename std::vector<std::string>::iterator itr =
		      gametes->second.begin() ; 
		    itr < gametes->second.end() ;
		    ++itr)
		  {
		    itr->resize(MAX_SEGSITES);
		  }
	      }
	    add_S_inf_sites(uni,i->begin(),tt,beg,end,i->nsam,nsites,
			    S,snp_index_i,gametes);
	    snp_index_i+=S;
	  }
      }
    std::sort(gametes->first.begin(),
	      gametes->first.begin()+ttlS);
    return ttlS;
  }

  template<typename uniform_generator>
  int infinite_sites_details( uniform_generator & uni,
			      gamete_storage_type * gametes,
			      const int & nsites,
			      const arg & history,
			      const double * total_times,
			      const unsigned * segsites )  
  {
    unsigned seg=0,nsegs = history.size();
    arg::const_iterator i=history.begin(),
      j = history.begin();
    j++;
    int ttlS = 0;
    for( ; seg < nsegs ; ++seg,++i,++j )
      {
	if( *(segsites+seg) > 0 )
	  {
	    bool add = false;
	    while(ttlS + int(*(segsites+seg)) > MAX_SEGSITES)
	      {
		MAX_SEGSITES += MAX_SEGS_INC;
		add=true;
	      }
	    if(add)
	      {
		gametes->first.resize(MAX_SEGSITES);
		for(typename std::vector<std::string>::iterator itr =
		      gametes->second.begin() ; 
		    itr < gametes->second.end() ;
		    ++itr)
		  {
		    itr->resize(MAX_SEGSITES);
		  }
	      }
	    int end = (seg<nsegs-1) ? j->beg : nsites;
	    int beg = i->beg;
	    add_S_inf_sites(uni,i->begin(),*(total_times+seg),
			    beg,end,i->nsam,nsites,
			    *(segsites+seg),ttlS,gametes);
	  }
	ttlS += *(segsites+seg);
      }
    std::sort(gametes->first.begin(),
	      gametes->first.begin()+ttlS);
    return ttlS;
  }



  template<typename poisson_generator,
	   typename uniform_generator>
  SimData infinite_sites_sim_data_details( poisson_generator & poiss,
					   uniform_generator & uni,
					   const int & nsites,
					   const arg & history,
					   const double & theta)
  {
    std::vector<unsigned> segsites(history.size());
    std::vector<double> total_times(history.size());
    arg::const_iterator i = history.begin(),j=i;
    j++;
    size_t nsegs = history.size();
    for(size_t seg = 0 ; seg < nsegs ; ++seg,++i,++j)
      {
	int end = (seg<nsegs-1) ? j->beg : nsites;
	int beg = i->beg;
	total_times[seg] = total_time(i->begin(),i->nsam);
	segsites[seg] = unsigned(poiss(boost::cref(total_times[seg]*theta*double(end-beg)/double(nsites))));
      }
    return infinite_sites_sim_data(uni,nsites,history,&total_times[0],&segsites[0]);
  }

  template<typename uniform_generator>
  SimData infinite_sites_sim_data_details( uniform_generator & uni,
					   const int & nsites,
					   const arg & history,
					   const double * total_times,
					   const unsigned * segsites )
  {
    assert(total_times != NULL);
    assert(segsites != NULL);
    const unsigned S = std::accumulate(segsites,segsites+history.size(),0u);
    if(S==0)
      return SimData();
    int nsam = history.begin()->nsam;
    SimData d(nsam,S);
    SimData::pos_iterator pos_i = d.pbegin();
    unsigned seg=0,nsegs = history.size();
    arg::const_iterator i=history.begin(),
      j = history.begin();
    j++;
    int ttlS = 0;
    //const double zero=0.;
    for( ; seg < nsegs ; ++seg,++i,++j )
      {
	if( *(segsites+seg) > 0 )
	  {
	    int end = (seg<nsegs-1) ? j->beg : nsites;
	    int beg = i->beg;
	    const double tt = *(total_times+seg);
	    for(unsigned snp = ttlS ; snp < ttlS+*(segsites+seg) ; ++snp)
	      {
		double pos = uni(beg,end)/double(nsites);
		*(pos_i+snp)=pos;
		int branch = pick_branch(i->begin(),nsam,uni(boost::cref(0.),tt));
		for(int ind=0;ind<nsam;++ind)
		  {
		    if( is_descendant(i->begin(),ind,branch) )
		      {
			d[ind][snp] = '1';
		      }
		    else
		      {
			d[ind][snp] = '0';
		      }
		  }
	      }
	    ttlS += *(segsites+seg);
	  }
      }
    std::sort(pos_i,d.pend());
    return d;
  }
#endif

  template<typename uniform_generator>
  void add_S_inf_sites ( uniform_generator & uni ,
			 marginal::const_iterator history,
			 const double & tt,
			 const int & beg, const int & end,
			 const int & nsam,
			 const int & nsites,
			 const int & S ,
			 const int & first_snp_index,
			 gamete_storage_type * gametes )
  /*! 
    @brief Add S segregating sites to sample with a particular marginal history, according
    to the infinitely-many sites model.  

    This routine places a fixed number of segregating  sites on a marginal history that
    begins at position \a beg and ends at position \a end-1.  
    Mutations are assigned positions randomly on the continuous, half-open interval 
    [ beg/L,end/L ), where L = nsites-1.

    \param uni a uniform random number generator that takes two doubles as arguments
    \param history the history onto which mutations will be placed
    \param tt the total time on \a history
    \param beg the first site in the history (0 <= beg < nsites)
    \param end the last site in the history ( beg < end < nsites )
    \param nsam the total sample size being simulated
    \param nsites the number of mutational sites simulated
    \param S the number of mutations to drop on \a history
    \param first_snp_index the index at which mutations will begin to be added into \a gametes
    \param gametes a container in which you are storing the mutations.  
    Must be allocated in the calling environment.
    \ingroup coalescent
  */
  {
    return add_S_inf_sites_details(uni,history,tt,beg,end,nsam,nsites,S,first_snp_index,gametes);
  }

  template<typename uniform_generator>
  void add_S_inf_sites ( const uniform_generator & uni ,
			 marginal::const_iterator history,
			 const double & tt,
			 const int & beg, const int & end,
			 const int & nsam,
			 const int & nsites,
			 const int & S ,
			 const int & first_snp_index,
			 gamete_storage_type * gametes )
  /*! 
    @brief Add S segregating sites to sample with a particular marginal history, according
    to the infinitely-many sites model.  

    This routine places a fixed number of segregating  sites on a marginal history that
    begins at position \a beg and ends at position \a end-1.  
    Mutations are assigned positions randomly on the continuous, half-open interval 
    [ beg/L,end/L ), where L = nsites-1.

    \param uni a uniform random number generator that takes two doubles as arguments
    \param history the history onto which mutations will be placed
    \param tt the total time on \a history
    \param beg the first site in the history (0 <= beg < nsites)
    \param end the last site in the history ( beg < end < nsites )
    \param nsam the total sample size being simulated
    \param nsites the number of mutational sites simulated
    \param S the number of mutations to drop on \a history
    \param first_snp_index the index at which mutations will begin to be added into \a gametes
    \param gametes a container in which you are storing the mutations.  
    Must be allocated in the calling environment.
    \ingroup coalescent
  */
  {
    return add_S_inf_sites_details(uni,history,tt,beg,end,nsam,nsites,S,first_snp_index,gametes);
  }

  template<typename poisson_generator,
	   typename uniform_generator>
  int infinite_sites( poisson_generator & poiss,
		      uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double & theta )
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph
    \param poiss a Poisson random number generator which takes the mean of the poisson as an argument
    \param uni a uniform random number generator that takes two doubles as an argument
    \param gametes object in which to store the simulated gametes
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param theta the coalescent-scaled mutation rate
    \return the number of mutations placed on the tree
    \ingroup coalescent
  */
  {
    return infinite_sites_details(poiss,uni,gametes,nsites,history,theta);
  }
  
  template<typename poisson_generator,
	   typename uniform_generator>
  int infinite_sites( const poisson_generator & poiss,
		      const uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double & theta )
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph
    \param poiss a Poisson random number generator which takes the mean of the poisson as an argument
    \param uni a uniform random number generator that takes two doubles as an argument
    \param gametes object in which to store the simulated gametes
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param theta the coalescent-scaled mutation rate
    \return the number of mutations placed on the tree
    \ingroup coalescent
  */     
  {
    return infinite_sites_details(poiss,uni,gametes,nsites,history,theta);
  }
  
  template<typename uniform_generator>
  int infinite_sites( uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double * total_times,
		      const unsigned * segsites )
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph with
    a fixed number of segregating sites
    \param uni a uniform random number generator that takes two doubles as an argument
    \param gametes object in which to store the simulated gametes
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param total_times the total times on each marginal tree in \a history
    \param segsites the number of segregating sites to place on each tree
    \return the number of mutations placed on the tree
    \note \a total_times and \a segsites must contain a number of elements equal to history.size()
    \ingroup coalescent
  */  
  {
    return infinite_sites_details(uni,gametes,nsites,history,total_times,segsites);
  } 

  template<typename uniform_generator>
  int infinite_sites( const uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double * total_times,
		      const unsigned * segsites )
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph with
    a fixed number of segregating sites
    \param uni a uniform random number generator that takes two doubles as an argument
    \param gametes object in which to store the simulated gametes
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param total_times the total times on each marginal tree in \a history
    \param segsites the number of segregating sites to place on each tree
    \return the number of mutations placed on the tree
    \note \a total_times and \a segsites must contain a number of elements equal to history.size()
    \ingroup coalescent
  */  

  {  return infinite_sites_details(uni,gametes,nsites,history,total_times,segsites);
  }

  template<typename poisson_generator,
	   typename uniform_generator>
  SimData infinite_sites_sim_data( poisson_generator & poiss,
				   uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double & theta)
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph
    \param poiss a Poisson random number generator which takes the mean of the poisson as an argument
    \param uni a uniform random number generator that takes two doubles as an argument
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param theta the coalescent-scaled mutation rate
    \return an object of type SimData that represent the sample. (the gametes are also
    stored in \a gametes).  The SimData object can be passed directly into class PolySIM
    for analysis
    \ingroup coalescent
  */
  {
    return infinite_sites_sim_data_details(poiss,uni,nsites,history,theta);
  }

  template<typename poisson_generator,
	   typename uniform_generator>
  SimData infinite_sites_sim_data( const poisson_generator & poiss,
				   const uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double & theta)
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph
    \param poiss a Poisson random number generator which takes the mean of the poisson as an argument
    \param uni a uniform random number generator that takes two doubles as an argument
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param theta the coalescent-scaled mutation rate
    \return an object of type SimData that represent the sample. (the gametes are also
    stored in \a gametes).  The SimData object can be passed directly into class PolySIM
    for analysis
    \ingroup coalescent
  */   
  {
    return infinite_sites_sim_data_details(poiss,uni,nsites,history,theta);
  }
 
  template<typename uniform_generator>
  SimData infinite_sites_sim_data( uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double * total_times,
				   const unsigned * segsites)
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph with
    a fixed number of segregating sites
    \param uni a uniform random number generator that takes two doubles as an argument
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param total_times the total times on each marginal tree in \a history
    \param segsites the number of segregating sites to place on each tree
    \return an object of type SimData that represent the sample. 
    for analysis
    \note \a total_times and \a segsites must contain a number of elements equal to history.size()
    \ingroup coalescent
  */
  {
    return infinite_sites_sim_data_details(uni,nsites,history,total_times,segsites);
  }

  template<typename uniform_generator>
  SimData infinite_sites_sim_data( const uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double * total_times,
				   const unsigned * segsites)
  /*!
    @brief Apply the infinitely-many sites mutation model to an ancetral recombination graph with
    a fixed number of segregating sites
    \param uni a uniform random number generator that takes two doubles as an argument
    \param nsites the length of the region begin simulated
    \param history the list of marginal histories for the sample
    \param total_times the total times on each marginal tree in \a history
    \param segsites the number of segregating sites to place on each tree
    \return an object of type SimData that represent the sample. 
    for analysis
    \note \a total_times and \a segsites must contain a number of elements equal to history.size()
    \ingroup coalescent
  */
  {
    return infinite_sites_sim_data_details(uni,nsites,history,total_times,segsites);
  }

} //ns sequence
#endif
