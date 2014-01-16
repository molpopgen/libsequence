#ifndef __SEQUENCE__COALESCENT__TRAJECTORIES_HPP__
#define __SEQUENCE__COALESCENT__TRAJECTORIES_HPP__

#include <vector>

namespace Sequence
{
  template< typename uni01_generator >
  void ConditionalTraj(uni01_generator & uni01,
		  std::vector<double> * traj,
		  const unsigned & N,
		  const double & s,
		  const double & dt,
		  const double & initial_frequency,
		  const double & final_frequency = 1.);
  
  template< typename uni01_generator >
  void ConditionalTraj(const uni01_generator & uni01,
		       std::vector<double> * traj,
		       const unsigned & N,
		       const double & s,
		       const double & dt,
		       const double & initial_frequency,
		       const double & final_frequency = 1.);

  template<typename uni01_generator>
  void ConditionalTrajNeutral(uni01_generator & uni01,
			      std::vector<double> * traj,
			      const double & dt,
			      const double & initial_freq = 1.,
			      const double & final_freq = 0.);
  
  template<typename uni01_generator>
  void ConditionalTrajNeutral(const uni01_generator & uni01,
			      std::vector<double> * traj,
			      const double & dt,
			      const double & initial_freq = 1.,
			      const double & final_freq = 0.);
}
#include <Sequence/Coalescent/bits/Trajectories.tcc>
#endif
