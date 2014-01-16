#ifndef __SEQUENCE_COALESCENT_DEMOGRAPHICMODELS_HPP__
#define __SEQUENCE_COALESCENT_DEMOGRAPHICMODELS_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>

namespace Sequence
{
  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg bottleneck( uniform_generator & uni,
		  uniform01_generator & uni01,
		  exponential_generator & expo,
		  const std::vector<chromosome> & initialized_sample, 
		  const marginal & initialized_marginal,
		  const double & tr,
		  const double & d,
		  const double & f,
		  const double & rho = 0.,
		  const bool & exponential_recovery = false,
		  const double & recovered_size = 1. );

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg bottleneck( const uniform_generator & uni,
		  const uniform01_generator & uni01,
		  const exponential_generator & expo,
		  const std::vector<chromosome> & initialized_sample, 
		  const marginal & initialized_marginal,
		  const double & tr,
		  const double & d,
		  const double & f,
		  const double & rho = 0.,
		  const bool & exponential_recovery = false,
		  const double & recovered_size = 1. );

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg exponential_change( uniform_generator & uni,
			  uniform01_generator & uni01,
			  exponential_generator & expo,
			  const std::vector<chromosome> & initialized_sample, 
			  const marginal & initialized_marginal,
			  const double & G,
			  const double & t_begin,
			  const double & t_end,
			  const double & rho = 0.,
			  const double & size_at_end = -1);

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg exponential_change( const uniform_generator & uni,
			  const uniform01_generator & uni01,
			  const exponential_generator & expo,
			  const std::vector<chromosome> & initialized_sample, 
			  const marginal & initialized_marginal,
			  const double & G,
			  const double & t_begin,
			  const double & t_end,
			  const double & rho = 0.,
			  const double & size_at_end = -1);

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg snm( uniform_generator & uni,
	   uniform01_generator & uni01,
	   exponential_generator & expo,
	   const std::vector<chromosome> & initialized_sample, 
	   const marginal & initialized_marginal,
	   const double & rho);

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg snm( const uniform_generator & uni,
	   const uniform01_generator & uni01,
	   const exponential_generator & expo,
	   const std::vector<chromosome> & initialized_sample, 
	   const marginal & initialized_marginal,
	   const double & rho);

} //namespace Sequence
#endif
#include <Sequence/Coalescent/bits/DemographicModels.tcc>
