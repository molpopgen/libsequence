// -*- C++ -*-
#ifndef __SEQUENCE_COALESCENT_BITS_TRAJECTORIES_TCC__
#define __SEQUENCE_COALESCENT_BITS_TRAJECTORIES_TCC__
#include <cassert>
#include <cmath>
#include <algorithm>
#include <boost/static_assert.hpp>

namespace Sequence
{
  template< typename uni01_generator >
  void ConditionalTraj_details(uni01_generator & uni01,
			       std::vector<double> * traj,
			       const unsigned & N,
			       const double & s,
			       const double & dt,
			       const double & initial_frequency,
			       const double & final_frequency)
  /*!
    Implementation details
  */
  {
    assert( initial_frequency > 0. );
    assert( final_frequency > 0. );
    assert( final_frequency <= 1. );
    assert( s>=0. );
    double x = initial_frequency;
    traj->erase(traj->begin(),traj->end());
    traj->push_back(x);
    while( x <= final_frequency )
      {
	double ux = (double(2*N)*s*x*(1.-x))/tanh(double(2*N)*s*x);
	if( uni01() < 0.5 )
	  {
	    x += ( ux*dt + std::sqrt( x*(1.-x)*dt ) );
	  }
	else
	  {
	    x += ( ux*dt - std::sqrt( x*(1.-x)*dt ) );
	  }
	assert( std::isfinite(x) );
	traj->push_back( ( (x<final_frequency)
			  ? x : final_frequency ) );
      }
    std::reverse(traj->begin(),traj->end());
  }

  template<typename uni01_generator>
  void ConditionalTrajNeutral_details(uni01_generator & uni01,
				      std::vector<double> * traj,
				      const double & dt,
				      const double & initial_freq,
				      const double & final_freq)
  /*!
    Implementation details
  */
  {
    double x = initial_freq;
    traj->erase(traj->begin(),traj->end());
    traj->push_back(x);
    while( x >= final_freq )
      {
	double ux = -1.*x;
	if( uni01() < 0.5 )
	  {
	    x += ( ux*dt + std::sqrt( x*(1.-x)*dt ) );
	  }
	else
	  {
	    x += ( ux*dt - std::sqrt( x*(1.-x)*dt ) );
	  }
	assert( std::isfinite(x) );
	traj->push_back( ( (x>final_freq) ? x : final_freq ) );
      }
  }

  template< typename uni01_generator >
  void ConditionalTraj(uni01_generator & uni01,
		       std::vector<double> * traj,
		       const unsigned & N,
		       const double & s,
		       const double & dt,
		       const double & initial_frequency,
		       const double & final_frequency)
  /*!
    Stochastic trajectory of beneficial mutations, following 
    Coop and Griffiths (2004).  
    \param uni01 a random number generator returning a U(0,1].
    \param traj For a diploid population of size N, this function will return trajectories
    equivalent to what one would get from a Wright-Fisher simulation of a haploid 
    population of size 2N
    \param dt amount by which to change increment time during simulation. 
    \param initial_frequency Initial frequency of beneficial allele (i.e. 1/2N)
    \param final_frequency Final frequency of beneficial allele (1 means fixation).
    \return A vector of length L, such that,for dt = 1/(k*2N), 
    L/(k*2N) is the length of the sweep, in units of 2N generations
  */
  {
    ConditionalTraj_details(uni01,traj,N,s,dt,initial_frequency,final_frequency);
  }
  
  template< typename uni01_generator >
  void ConditionalTraj(const uni01_generator & uni01,
		       std::vector<double> * traj,
		       const unsigned & N,
		       const double & s,
		       const double & dt,
		       const double & initial_frequency,
		       const double & final_frequency)
  /*!
    Stochastic trajectory of beneficial mutations, following 
    Coop and Griffiths (2004).  
    \param uni01 a random number generator returning a U(0,1].
    \param traj A vector of length L, such that,for dt = 1/(k*2N), 
    L/(k*2N) is the length of the sweep, in units of 2N generations
    \note For a diploid population of size N, this function will return trajectories
    equivalent to what one would get from a Wright-Fisher simulation of a haploid 
    population of size 2N
    \param dt amount by which to change increment time during simulation. 
    \param initial_frequency Initial frequency of beneficial allele (i.e. 1/2N)
    \param final_frequency Final frequency of beneficial allele (1 means fixation).
  */
  {
    ConditionalTraj_details(uni01,traj,N,s,dt,initial_frequency,final_frequency);
  }
  
  template<typename uni01_generator>
  void ConditionalTrajNeutral(uni01_generator & uni01,
			      std::vector<double> * traj,
			      const double & dt,
			      const double & initial_freq,
			      const double & final_freq)
  
  /*!
    Stochastic trajectory of a neutral allele, following Coop & Griffiths (2004) TPB,
    and Przeworski et al. (2005) Evolution.  The simulation is backwards in time.
    \param uni01 a random number generator returning a U(0,1].
    \param traj A vector of length L, and for dt=1/(k*2N), the length of time,
    in units of 2N generations, is given by L/(k*2N).  The vector describes
    the change in allele frequency from \a initial_freq to \a final_freq,
    in jumps in time of \a dt.
    \param dt amount by which to change increment time during simulation. 
    \param initial_freq Initial frequency of the neutral allele.
    \param final_freq final frequency of the neutral allele.
  */
  {
    ConditionalTrajNeutral_details(uni01,traj,dt,initial_freq,final_freq);
  }
  
  template<typename uni01_generator>
  void ConditionalTrajNeutral(const uni01_generator & uni01,
			      std::vector<double> * traj,
			      const double & dt,
			      const double & initial_freq,
			      const double & final_freq)
  /*!
    Stochastic trajectory of a neutral allele, following Coop & Griffiths (2004) TPB,
    and Przeworski et al. (2005) Evolution.  The simulation is backwards in time.
    \param uni01 a random number generator returning a U(0,1].
    \param traj A vector of length L, and for dt=1/(k*2N), the length of time,
    in units of 2N generations, is given by L/(k*2N).  The vector describes
    the change in allele frequency from \a initial_freq to \a final_freq,
    in jumps in time of \a dt.
    \param dt amount by which to change increment time during simulation. 
    \param initial_freq Initial frequency of the neutral allele.
    \param final_freq final frequency of the neutral allele.
  */
  {
    ConditionalTrajNeutral_details(uni01,traj,dt,initial_freq,final_freq);
  }
}
#endif
