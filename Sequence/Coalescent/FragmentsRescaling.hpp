/*! \file FragmentsRescaling.hpp
  \brief Helper functions for simulating partially linked fragments
  One often wants to simulate partially linked fragments under neutral models.
  An efficient way to do this is to simulate a contiguous fragment, but
  with the recombination rate varying along the region (to represent
  the variable genetic distances between fragments).  This header file
  declares functions that make this task easier, particularly the 
  operations of rescaling the positions of mutations/marginal trees
  from the genetic map back to the physical map
*/

#ifndef __SEQUENCE_COALESCENT_FRAGMENTS_RESCALING_HPP__
#define __SEQUENCE_COALESCENT_FRAGMENTS_RESCALING_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <vector>
#include <utility>

namespace Sequence
{
  int sample_length( const std::vector< std::pair<int,int> > & fragments );
  int total_length( const std::vector< std::pair<int,int> > & fragments );
  void calculate_scales(const std::vector< std::pair<int,int> > & fragments,
			std::vector< std::pair<double,double> > * sample_scale,
			std::vector< std::pair<double,double> > * mutation_scale );
  class SimData; //forward declaration
  void rescale_mutation_positions(SimData * d,
				  const std::vector< std::pair<double,double> > & sample_scale, 
				  const std::vector< std::pair<double,double> > & mutation_scale );
  void rescale_arg( arg * sample_history,
		    const std::vector< std::pair<int,int> > & fragments );
  double integrate_genetic_map( const std::vector<chromosome> & sample,
				const int & current_nsam,
				const std::vector<double> & genetic_map,
				std::vector<double> * reclens);
}

#endif
