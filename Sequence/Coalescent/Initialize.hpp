#ifndef __SEQUENCE_COALESCENT_INIT_ARG_FUNCTIONS_HPP_
#define __SEQUENCE_COALESCENT_INIT_ARG_FUNCTIONS_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <vector>

namespace Sequence
{
  std::vector<chromosome> init_sample( const std::vector<int> & pop_config,
				       const int & nsites );
  marginal init_marginal( const int & nsam );
}

#endif
