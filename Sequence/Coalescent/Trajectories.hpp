#ifndef __SEQUENCE__COALESCENT__TRAJECTORIES_HPP__
#define __SEQUENCE__COALESCENT__TRAJECTORIES_HPP__

#include <vector>

namespace Sequence
{
  namespace coalsim {
    /*!
      Generate a trajectory for an additive beneficial mutation
      with selection coefficient "s", follwing the method of 
      Coop and Griffiths (http://www.ncbi.nlm.nih.gov/pubmed/15465123)

      @param uni01 Generates doubles on the interval [0,1)
      @param traj This will contain the trjactory.  It is cleared by the function.
      @param N The population size (keep it small, <= 10^4 or so)
      @param s The selection coefficient
      @param dt Determines the time scale
      @param initial_frequency The starting frequency of the beneficial mutation
      @param final_frequency The final_frequency of the beneficial mutation

      \note Upon return, traj contains frequency values starting at @a final_frequency and ending
      at @a initial_frequency
     */
    template< typename uni01_generator >
    void ConditionalTraj(uni01_generator & uni01,
			 std::vector<double> * traj,
			 const unsigned & N,
			 const double & s,
			 const double & dt,
			 const double & initial_frequency,
			 const double & final_frequency = 1.);
  
    /*!
      Generate a trajectory for an additive beneficial mutation
      with selection coefficient "s", follwing the method of 
      Coop and Griffiths (http://www.ncbi.nlm.nih.gov/pubmed/15465123)

      @param uni01 Generates doubles on the interval [0,1)
      @param traj This will contain the trjactory.  It is cleared by the function.
      @param N The population size (keep it small, <= 10^4 or so)
      @param s The selection coefficient
      @param dt Determines the time scale
      @param initial_frequency The starting frequency of the beneficial mutation
      @param final_frequency The final_frequency of the beneficial mutation

      \note Upon return, traj contains frequency values starting at @a final_frequency and ending
      at @a initial_frequency
     */
    template< typename uni01_generator >
    void ConditionalTraj(const uni01_generator & uni01,
			 std::vector<double> * traj,
			 const unsigned & N,
			 const double & s,
			 const double & dt,
			 const double & initial_frequency,
			 const double & final_frequency = 1.);

    /*!
      Generate a trajectory for a neutral mutation via the method
      of Przeworski et al. (2005),
      http://dx.doi.org/10.1554/05-273.1

      @param uni01 Generates doubles on the interval [0,1)
      @param traj This will contain the trjactory.  It is cleared by the function.
      @param dt Determines the time scale
      @param initial_frequency The starting frequency of the beneficial mutation
      @param final_frequency The final_frequency of the beneficial mutation

      \note Upon return, traj contains frequency values starting at @a final_frequency and ending
      at @a initial_frequency
    */
    template<typename uni01_generator>
    void ConditionalTrajNeutral(uni01_generator & uni01,
				std::vector<double> * traj,
				const double & dt,
				const double & initial_freq = 1.,
				const double & final_freq = 0.);
  
    /*!
      Generate a trajectory for a neutral mutation via the method
      of Przeworski et al. (2005),
      http://dx.doi.org/10.1554/05-273.1

      @param uni01 Generates doubles on the interval [0,1)
      @param traj This will contain the trjactory.  It is cleared by the function.
      @param dt Determines the time scale
      @param initial_frequency The starting frequency of the beneficial mutation
      @param final_frequency The final_frequency of the beneficial mutation

      \note Upon return, traj contains frequency values starting at @a final_frequency and ending
      at @a initial_frequency
    */
    template<typename uni01_generator>
    void ConditionalTrajNeutral(const uni01_generator & uni01,
				std::vector<double> * traj,
				const double & dt,
				const double & initial_freq = 1.,
				const double & final_freq = 0.);
  }
}
#include <Sequence/Coalescent/bits/Trajectories.tcc>
#endif
