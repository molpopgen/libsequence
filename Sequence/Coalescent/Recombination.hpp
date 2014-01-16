#ifndef __SEQUENCE_COALESCENT_RECOMBINATION_HPP__
#define __SEQUENCE_COALESCENT_RECOMBINATION_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
namespace Sequence
{
  int crossover( const int & current_nsam,
		 const int & chromo,
		 const int & pos,
		 std::vector<chromosome> * sample,
		 arg * sample_history);

  std::pair<int,int> pick_uniform_spot(const double & random_01,
				       const int & nlinks,
				       std::vector<chromosome>::const_iterator sample_begin,
				       const unsigned & current_nsam);

  template<typename uniform01_generator>
  std::pair<int,int> pick_spot( uniform01_generator & uni01,
				const double & total_reclen,
				const std::vector<double> & reclens,
				std::vector<chromosome>::const_iterator sample_begin,
				const unsigned & current_nsam,
				const double * rec_map);

  template<typename uniform01_generator>
  std::pair<int,int> pick_spot( const uniform01_generator & uni01,
				const double & total_reclen,
				const std::vector<double> & reclens,
				std::vector<chromosome>::const_iterator sample_begin,
				const unsigned & current_nsam,
				const double * rec_map);
}
#endif
#include <Sequence/Coalescent/bits/Recombination.tcc>
