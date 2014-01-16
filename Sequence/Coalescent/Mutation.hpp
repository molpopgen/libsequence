#ifndef __SEQUENCE_COALESCENT_MUTATION_HPP__
#define __SEQUENCE_COALESCENT_MUTATION_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/SimData.hpp>
#include <vector>
#include <string>
#include <utility>
#include <cstdio>

namespace Sequence
{
  /*! 
    @brief controls allocation of simulated gametes
    You must define this in namespace Sequence in your program.
    A value of 200 works well.
    \ingroup coalescent
  */
  extern int MAX_SEGSITES;
  /*! 
    @brief controls (re)allocation of simulated gametes
    You must define this in namespace Sequence in your program.
    A value of 100 works well
    \ingroup coalescent
  */
  extern int MAX_SEGS_INC;

  /*!
    @brief an object to store simulated gametes
    An object of this type will tend to exist in the calling environment
    of your program. If you are simulating a sample of n chromosomes, you
    would initialize the object as follows:
    \code
    gamete_storage_type gamete_bucket( std::vector<double>(MAX_SEGSITES,0.),
    std::vector< std::string >(n,std::string(MAX_SEGSITES,'0')) );
    \endcode
    \ingroup coalescent
   */
  typedef std::pair< std::vector<double>, 
		     std::vector<std::string> > gamete_storage_type;

  template<typename uniform_generator>
  void add_S_inf_sites ( uniform_generator & uni ,
			 marginal::const_iterator history,
			 const double & tt,
			 const int & beg, const int & end,
			 const int & nsam,
			 const int & nsites,
			 const int & S ,
			 const int & first_snp_index,
			 gamete_storage_type * gametes );

  template<typename uniform_generator>
  void add_S_inf_sites ( const uniform_generator & uni,
			 marginal::const_iterator history,
			 const double & tt,
			 const int & beg, const int & end,
			 const int & nsam,
			 const int & nsites,
			 const int & S ,
			 const int & first_snp_index,
			 gamete_storage_type * gametes );

  template<typename poisson_generator,
	   typename uniform_generator>
  int infinite_sites( poisson_generator & poiss,
		      uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double & theta );

  template<typename poisson_generator,
	   typename uniform_generator>
  int infinite_sites( const poisson_generator & poiss,
		      const uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double & theta );

  template<typename uniform_generator>
  int infinite_sites( uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double * total_times,
		      const unsigned * segsites );  

  template<typename uniform_generator>
  int infinite_sites( const uniform_generator & uni,
		      gamete_storage_type * gametes,
		      const int & nsites,
		      const arg & history,
		      const double * total_times,
		      const unsigned * segsites );  

  template<typename poisson_generator,
	   typename uniform_generator>
  SimData infinite_sites_sim_data( poisson_generator & poiss,
				   uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double & theta);

  template<typename poisson_generator,
	   typename uniform_generator>
  SimData infinite_sites_sim_data( const poisson_generator & poiss,
				   const uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double & theta);
 
  template<typename uniform_generator>
  SimData infinite_sites_sim_data( uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double * total_times,
				   const unsigned * segsites);

  template<typename uniform_generator>
  SimData infinite_sites_sim_data( const uniform_generator & uni,
				   const int & nsites,
				   const arg & history,
				   const double * total_times,
				   const unsigned * segsites);

  void output_gametes(FILE * fp,const unsigned & segsites,
		      const unsigned & nsam,
		      const gamete_storage_type & gametes);
}
#endif
#include <Sequence/Coalescent/bits/Mutation.tcc>
