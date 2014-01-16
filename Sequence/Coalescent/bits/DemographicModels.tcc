//  -*- C++ -*- 
#ifndef __SEQUENCE_COALESCENT_BITS_DEMOGRAPHIC_MODELS_HPP__
#define __SEQUENCE_COALESCENT_BITS_DEMOGRAPHIC_MODELS_HPP__

#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <Sequence/SeqConstants.hpp>
#include <boost/ref.hpp>
#include <limits>
#include <cmath>
#include <cassert>

namespace Sequence
{
#ifndef DOXYGEN_SKIP
  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg bottleneck_details( uniform_generator & uni,
			  uniform01_generator & uni01,
			  exponential_generator & expo,
			  const std::vector<chromosome> & initialized_sample, 
			  const marginal & initialized_marginal,
			  const double & tr,
			  const double & d,
			  const double & f,
			  const double & rho,
			  const bool & exponential_recovery,
			  const double & recovered_size )
  {
    assert( d > 0. );
    assert( tr > 0. );
    assert( f > 0. );
    assert( rho >= 0. );
    assert( recovered_size > 0. );
    assert( initialized_marginal.nsam == int(initialized_sample.size()) );

    arg sample_history(1,initialized_marginal);
    std::vector<chromosome> sample(initialized_sample);
    int nsam=sample.size(),NSAM = sample.size();
    const int nsites = sample[0].last()+1;
    const double littler = rho/double(nsites-1);
    int nlinks = NSAM*(nsites-1);
    double t = 0.;
    const double G = (std::log(recovered_size)-std::log(f))/d;
    double tcoal,rcoal,trec,rrec,tmin;
    bool event;
    std::pair<int,int> two;
    while( NSAM > 1 )
      {
	rcoal = double(NSAM*(NSAM-1));
	rrec = littler*(nlinks);
	if ( exponential_recovery )
	  {
	    assert(rcoal>0.);
	    if(t<tr)
	      {
		tcoal = expo(boost::cref(recovered_size/rcoal));
	      }
	    if(t >= tr && t < (tr+d))
	      {
		//uni01() must return an rv ~ U[0,1), so we must use 1-rv for logs...
		double temp = 1.-G*recovered_size*std::exp(-G*(t-tr))*std::log(1.-uni01())/rcoal;
		assert(temp >= 0.);
		tcoal = std::log(temp)/G;
		assert(std::isfinite(tcoal));
	      }
	    else
	      {
		tcoal = expo(boost::cref(1./rcoal));
	      }
	  }
	else
	  {
	    assert(rcoal > 0.);
	    if(t<tr)
	      {
		tcoal = expo(boost::cref(recovered_size/rcoal));
	      }
	    if(t >= tr && t < (tr+d))
	      {
		tcoal = expo(boost::cref(f/rcoal));
	      }
	    else
	      {
		tcoal = expo(boost::cref( 1./rcoal));
	      }
	  }
	trec = (rrec>0.) ? expo(boost::cref(1./rrec)) : SEQMAXDOUBLE;
	assert(std::isfinite(trec));
	tmin = std::min(tcoal,trec);
	event = true;
	if( t < tr && (t+tmin) >= tr )
	  {
	    t = tr;
	    event=false;
	  }
	else if ( t < (tr+d) && (t+tmin) >= (tr+d) )
	  {
	    t = tr+d;
	    event=false;
	  }
	if(event)
	  {
	    if (tcoal < trec)
	      {
		t += tmin;
		two = pick2(uni,NSAM);
		NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
				 &nlinks,&sample,&sample_history);
	      }
	    else
	      {
		t += tmin;
		two = pick_uniform_spot(uni01(),
					nlinks,
					sample.begin(),NSAM);
		nlinks -= crossover(NSAM,two.first,two.second,
				    &sample,&sample_history);
		NSAM++;
	      }
	    if(NSAM < int(sample.size())/5)
	      sample.erase(sample.begin()+NSAM+1,sample.end());
	  }
      }
    return sample_history;
  }

  //works and tested against ms
  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg exponential_change_details( uniform_generator & uni,
				  uniform01_generator & uni01,
				  exponential_generator & expo,
				  const std::vector<chromosome> & initialized_sample, 
				  const marginal & initialized_marginal,
				  const double & G,
				  const double & t_begin,
				  const double & t_end,
				  const double & rho,
				  const double & size_at_end)
  {
    assert( t_begin >= 0. );
    assert( t_end >= 0. );
    assert( t_end >= t_begin );
    assert( rho >= 0. );
    assert( initialized_marginal.nsam == int(initialized_sample.size()) );
    arg sample_history(1,initialized_marginal);
    std::vector<chromosome> sample(initialized_sample);
    int nsam=sample.size(),NSAM = sample.size();
    const int nsites = sample[0].last()+1;
    const double littler = rho/double(nsites-1);
    int nlinks = NSAM*(nsites-1);
    double t = 0.,tcoal,trec,rrec,rcoal,tmin;
    bool event;
    std::pair<int,int> two;
    double temp;
    while( NSAM > 1 )
      {
	assert( std::isfinite(t) );
	rcoal = double(NSAM*(NSAM-1));
	rrec = littler*(nlinks);
	trec = (rrec>0.) ? expo(boost::cref(1./rrec)) : SEQMAXDOUBLE;
	assert(std::isfinite(trec));
	if ( t < t_begin || G == 0.)
	  {
	    tcoal = expo(boost::cref(1./rcoal));
	    assert( std::isfinite(tcoal) );
	  }
	else if ( t >= t_begin && t < t_end )
	  {
	    temp = 1.-G*std::exp(-G*(t-t_begin))*std::log(1.-uni01())/rcoal;
	    assert( temp >= 0. );
	    tcoal = std::log(temp)/G;
	    assert( std::isfinite(tcoal) );
	  }
	else
	  {
	    tcoal = (size_at_end > 0.) ? expo(boost::cref(size_at_end/rcoal)): 
	      expo(boost::cref(std::exp(-G*(t_end-t_begin))/rcoal));
	    assert( std::isfinite(tcoal) );
	  }
	tmin = std::min(tcoal,trec);
	assert( std::isfinite(tmin) );
	event = true;
	if( t < t_begin && (t+tmin) >= t_begin)
	  {
	    t=t_begin;
	    event=false;
	  }
	else if ( t<t_end && (t+tmin) >= t_end )
	  {
	    t=t_end;
	    event=false;
	  }
	if(event)
	  {
	    if (tcoal < trec)
	      {
		t += tmin;
		two = pick2(uni,NSAM);
		NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
				 &nlinks,&sample,&sample_history);
	      }
	    else
	      {
		t += tmin;
		two = pick_uniform_spot(uni01(),
					nlinks,
					sample.begin(),NSAM);
		nlinks -= crossover(NSAM,two.first,two.second,
				    &sample,&sample_history);
		NSAM++;
	      }
	  }
	if(NSAM < int(sample.size())/5)
	  sample.erase(sample.begin()+NSAM+1,sample.end());
      }
    return sample_history; 
  }
#endif

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
		  const double & rho,
		  const bool & exponential_recovery,
		  const double & recovered_size )
  /*!
    @brief Coalescent simulation of a population bottleneck
    Simulate a single, bottlenecked, population according to the Wright-Fisher model without selection. The population can recover from the bottleneck either instantaneously ("stepwise bottleneck"), or according to an exponential growth model. For the case of a stepwise bottleneck, this function is equivalent to the following options in Dick Hudson's program "ms": -eN 0 recovered_size -eN tr f -eN (tr+d) 1.  For the case where recovery from the bottleneck is by exponential growth, the equivalent "ms" options are: -eN 0 recovered_size -eG tr (log(recovered_size)-log(f))/d -eG (tr+d) 0 -eN (tr+d) 1.
    \param uni A binary function object (or equivalent) that returns a random deviate between a and b such that a <= x < b.  a and b are the arguments to operator() of \a uni
    \param uni01 A function object (or equivalent) whose operator() takes no arguments and returns a random deviate 0 <= x < 1.
    \param expo A unary function object whose operator() takes the mean of an exponential process as an argument and returns a deviate from an exponential distribution with that mean
    \param initialized_sample An initialized vector of chromosomes for a single population.  For example, this may be the return value of init_sample.  This object is used to copy-construct a non-const sample for the simulation
    \param initialized_marginal  An initialized marginal tree of the appropriate sample size for the simulation.  For example, the return value of init_marginal.
    \param tr The time at which the population recovers from the bottleneck. In units of 4N0 generations, where N0 is the effective size before the bottleneck.
    \param d The duration of the bottleneck, in units of 4N0 generations, where N0 is the effective size before the bottleneck.
    \param f Bottleneck severity.  Define Nb as the effective size during the bottleneck, and N0 the effective size prior to the bottleneck.  f=Nb/N0.
    \param rho The population recombination rate 4N0r.  The number of "sites" simulated is not neccesary, as it can be obtained from initialized_sample[0].last()+1.
    \param exponential_recovery If true, the population recovers from the bottleneck according to an exponential growth model.  If false, a stepwise bottleneck is assumed.
    \param recovered_size  If 1, the population recovers to N0 at time \a tr.  If 0.5, the population recovers to 1/2 the pre-bottleneck size, etc.
    \return The ancestral recombination graph (arg) describing the sample history.
    \pre d>0 and tr>0 and f>0 and rho>=0 and recovered_size>0 and initialized_marginal.nsam == initialized_sample.size()
    \note Preconditions are checked by the assert macro, and are therefore disabled when compiling with -DNDEBUG.  No checks or warnings are otherwised performed nor given.
    \ingroup coalescent
  */
  {
    return bottleneck_details(uni,uni01,expo,initialized_sample,initialized_marginal,tr,d,f,rho,exponential_recovery,recovered_size);
  }

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
		  const double & rho,
		  const bool & exponential_recovery,
		  const double & recovered_size)
  /*!
    @brief Coalescent simulation of a population bottleneck
    Simulate a single, bottlenecked, population according to the Wright-Fisher model without selection. The population can recover from the bottleneck either instantaneously ("stepwise bottleneck"), or according to an exponential growth model. For the case of a stepwise bottleneck, this function is equivalent to the following options in Dick Hudson's program "ms": -eN 0 recovered_size -eN tr f -eN (tr+d) 1.  For the case where recovery from the bottleneck is by exponential growth, the equivalent "ms" options are: -eN 0 recovered_size -eG tr (log(recovered_size)-log(f))/d -eG (tr+d) 0 -eN (tr+d) 1.
    \param uni A binary function object (or equivalent) that returns a random deviate between a and b such that a <= x < b.  a and b are the arguments to operator() of \a uni
    \param uni01 A function object (or equivalent) whose operator() takes no arguments and returns a random deviate 0 <= x < 1.
    \param expo A unary function object whose operator() takes the mean of an exponential process as an argument and returns a deviate from an exponential distribution with that mean
    \param initialized_sample An initialized vector of chromosomes for a single population.  For example, this may be the return value of init_sample.  This object is used to copy-construct a non-const sample for the simulation
    \param initialized_marginal  An initialized marginal tree of the appropriate sample size for the simulation.  For example, the return value of init_marginal.
    \param tr The time at which the population recovers from the bottleneck. In units of 4N0 generations, where N0 is the effective size before the bottleneck.
    \param d The duration of the bottleneck, in units of 4N0 generations, where N0 is the effective size before the bottleneck.
    \param f Bottleneck severity.  Define Nb as the effective size during the bottleneck, and N0 the effective size prior to the bottleneck.  f=Nb/N0.
    \param rho The population recombination rate 4N0r.  The number of "sites" simulated is not neccesary, as it can be obtained from initialized_sample[0].last()+1.
    \param exponential_recovery If true, the population recovers from the bottleneck according to an exponential growth model.  If false, a stepwise bottleneck is assumed.
    \param recovered_size  If 1, the population recovers to N0 at time \a tr.  If 0.5, the population recovers to 1/2 the pre-bottleneck size, etc.
    \return The ancestral recombination graph (arg) describing the sample history.
    \pre d>0 and tr>0 and f>0 and rho>=0 and recovered_size>0 and initialized_marginal.nsam == initialized_sample.size()
    \note Preconditions are checked by the assert macro, and are therefore disabled when compiling with -DNDEBUG.  No checks or warnings are otherwised performed nor given.
    \ingroup coalescent
  */
  {
    return bottleneck_details(uni,uni01,expo,initialized_sample,initialized_marginal,tr,d,f,rho,exponential_recovery,recovered_size);
  }

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
			  const double & rho,
			  const double & size_at_end)
  /*!
    @brief Coalescent simulation of exponential change in population size
    Simulate a single population whose size changes exponentially during some period of time.  The relevant command line options for Hudson's program "ms" would be: -eG t_begin G -eG t_end 0. -eN t_end size_at_end.
    \param uni A binary function object (or equivalent) that returns a random deviate between a and b such that a <= x < b.  a and b are the arguments to operator() of \a uni
    \param uni01 A function object (or equivalent) whose operator() takes no arguments and returns a random deviate 0 <= x < 1.
    \param expo A unary function object whose operator() takes the mean of an exponential process as an argument and returns a deviate from an exponential distribution with that mean
    \param initialized_sample An initialized vector of chromosomes for a single population.  For example, this may be the return value of init_sample.  This object is used to copy-construct a non-const sample for the simulation
    \param initialized_marginal  An initialized marginal tree of the appropriate sample size for the simulation.  For example, the return value of init_marginal.
    \param G The rate of exponential change in effective size.  If G>0, the population grows exponentially (forwards in time).  If G<0, it shrinks (again, forwards in time).
    \param t_begin The time in the past (in units of 4Ne generations) at which population size change begins (i.e., ends, moving forward in time)
    \param t_end The time in the past (in units of 4Ne generations) at which populations size change ends (begins forward in time)
    \param rho The population recombination rate 4N0r.  The number of "sites" simulated is not neccesary, as it can be obtained from initialized_sample[0].last()+1.
    \param size_at_end  At time \a t_end in the past, the population size is set to \a size_at_end.  If \a size_at_end = 1, the population is set to the same size that is was at t=0 (i.e. the beginning of the simulation).  If \a size_at_and < 0, the population size is not adjusted at \a t_end.  In other words, it is left at whatever it grew or shrank to.
    \return The ancestral recombination graph (arg) describing the sample history.
    \pre t_begin>=0 and t_end>=0 and t_end>=t_begin and rho>=0  and initialized_marginal.nsam == initialized_sample.size()
    \note Preconditions are checked by the assert macro, and are therefore disabled when compiling with -DNDEBUG.  No checks or warnings are otherwised performed nor given.
    \ingroup coalescent
  */
  {
    return exponential_change_details(uni,uni01,expo,initialized_sample,initialized_marginal,
				      G,t_begin,t_end,rho,size_at_end);
  }

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
			  const double & rho,
			  const double & size_at_end)
  /*!
    @brief Coalescent simulation of exponential change in population size
    Simulate a single population whose size changes exponentially during some period of time.  The relevant command line options for Hudson's program "ms" would be: -eG t_begin -eG t_end 0. -eN t_end size_at_end.
    \param uni A binary function object (or equivalent) that returns a random deviate between a and b such that a <= x < b.  a and b are the arguments to operator() of \a uni
    \param uni01 A function object (or equivalent) whose operator() takes no arguments and returns a random deviate 0 <= x < 1.
    \param expo A unary function object whose operator() takes the mean of an exponential process as an argument and returns a deviate from an exponential distribution with that mean
    \param initialized_sample An initialized vector of chromosomes for a single population.  For example, this may be the return value of init_sample.  This object is used to copy-construct a non-const sample for the simulation
    \param initialized_marginal  An initialized marginal tree of the appropriate sample size for the simulation.  For example, the return value of init_marginal.
    \param G The rate of exponential change in effective size.  If G>0, the population grows exponentially (forwards in time).  If G<0, it shrinks (again, forwards in time).
    \param t_begin The time in the past (in units of 4Ne generations) at which population size change begins (i.e., ends, moving forward in time)
    \param t_end The time in the past (in units of 4Ne generations) at which populations size change ends (begins forward in time)
    \param rho The population recombination rate 4N0r.  The number of "sites" simulated is not neccesary, as it can be obtained from initialized_sample[0].last()+1.
    \param size_at_end  At time \a t_end in the past, the population size is set to \a size_at_end.  If \a size_at_end = 1, the population is set to the same size that is was at t=0 (i.e. the beginning of the simulation).  If \a size_at_and < 0, the population size is not adjusted at \a t_end.  In other words, it is left at whatever it grew or shrank to.
    \return The ancestral recombination graph (arg) describing the sample history.
    \pre t_begin>=0 and t_end>=0 and t_end>=t_begin and rho>=0  and initialized_marginal.nsam == initialized_sample.size()
    \note Preconditions are checked by the assert macro, and are therefore disabled when compiling with -DNDEBUG.  No checks or warnings are otherwised performed nor given.
    \ingroup coalescent
  */
  {
    return exponential_change_details(uni,uni01,expo,initialized_sample,initialized_marginal,
				      G,t_begin,t_end,rho,size_at_end);
  }

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg snm( uniform_generator & uni,
	   uniform01_generator & uni01,
	   exponential_generator & expo,
	   const std::vector<chromosome> & initialized_sample, 
	   const marginal & initialized_marginal,
	   const double & rho)
  {
    return exponential_change_details(uni,uni01,expo,initialized_sample,initialized_marginal,
				      0.,0.,0.,rho,1.);
  }

  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator>
  arg snm( const uniform_generator & uni,
	   const uniform01_generator & uni01,
	   const exponential_generator & expo,
	   const std::vector<chromosome> & initialized_sample, 
	   const marginal & initialized_marginal,
	   const double & rho)
  {
    return exponential_change_details(uni,uni01,expo,initialized_sample,initialized_marginal,
				      0.,0.,0.,rho,1.);
  }
} //namespace Sequence

#endif //include guard
