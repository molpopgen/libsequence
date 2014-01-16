#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <Sequence/Coalescent/Mutation.hpp>
#include <Sequence/SimData.hpp>
#include <utility>

namespace Sequence
{						
  template<typename uniform_generator,
	   typename uniform01_generator,
	   typename exponential_generator,
	   typename poisson_generator>
  Sequence::SimData neutral_sample( uniform_generator & uni,
				    uniform01_generator & uni01,
				    exponential_generator & expo,
				    poisson_generator & poiss,
				    const double & theta,
				    const double & rho,
				    const int & nsites,
				    const int & nsam,
				    std::vector<chromosome> * sample,
				    arg * sample_history,
				    unsigned * max_chromosomes = NULL,
				    const unsigned & max_chromosomes_inc = 0)
  /*!
    @brief A simple function to generate samples under a neutral equilibrium model.

    A simple function to generate samples under a neutral equilibrium model
    with infinite-sites mutation and a constant recombination rate accross the region.

    \param uni a function/object capable of returning a random double uniformly from [0,k)
    \param uni01 a function/object capable of returning a random probability uniformly from [0,1)
    \param expo a function/object capable of returning an exponentially distributed random variable. 
    The function must take a single double as an argument, which is the mean of the exponential 
    distribution
    \param poiss a function/object capable of returning an poisson distributed random variable. 
    The function must take a single double as an argument, which is the mean of the poisson 
    distribution
    \param theta 4Nu, the coalescent-scaled mutation rate
    \param rho 4Nr, the recombination rate for the whole region
    \param nsites the number of mutational sites to simulate.  Recombination is equally likely
    between any two sites.
    \param nsites the total sample size. (There is no population structure in this routine)
    \param sample A pointer to the sample of chromosomes you wish to simulate.  
    This must be properly initialized, for example using the function init_sample in 
    <Sequence/Coalescent/Initialize.hpp>
    \param sample_history a pointer to the ancestral recombination graph.  This must be
    initialized in the calling enviroment.  In general, you can use init_marginal in 
    <Sequence/Coalescent/Initialize.hpp>
    \param max_chromosomes  This is a pointer to an integer in the calling environment which you
    can use to reserve memory in the array containing the sample of chromosomes.  If the size of \a sample
    ever gets larger than this, max_chromosomes is incremented by \a max_chromosomes_inc
    \param max_chromosomes_inc the amount by which to increment \a max_chromosomes
    \note This function does require a bit of work to use, although not much.  Please see the example
    code that comes with the library, in particular ms--.cc
    \ingroup coalescent
   */
  {
    int NSAM = nsam;

    //this is rho = 4Nr/site
    double littler = rho/(double(nsites-1));

    //a chromosome with nsites sites has nsites-1 positions ("links")
    //at which crossovers can occur, so the total number of
    //links at the start of the simulation is:
    int nlinks = nsam*(nsites-1);
    double t = 0.;
    while(NSAM>1)
      {
	double rcoal = double(NSAM*(NSAM-1));
	double rrec = (rho>0.) ? littler*double(nlinks) : 0.;

	//note--the function calls below scale time in units of 4Ne
	double tcoal = expo(1./rcoal);
	double trec = expo(1./rrec);
	if ( trec < tcoal ) //crossover event
	  {
	    t+=trec;
	    std::pair<int,int> pos_rec = pick_uniform_spot(uni01(),
							   nlinks,
							   sample->begin(),NSAM);
 	    assert( pos_rec.second >= 0 );
 	    assert( pos_rec.second >= (sample->begin()+pos_rec.first)->first() );
	    assert( pos_rec.second <= ((sample->begin()+pos_rec.first)->last() ) ); 
				       
	    assert( (sample->begin()+pos_rec.first)->links()>0 );
	    nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
				sample,sample_history);
	    NSAM++;
	  }
	else //common ancestor event
	  {
	    t+=tcoal;
	    std::pair<int,int> two = pick2(uni,NSAM);
	    NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
			     &nlinks,sample,sample_history);
	  }
	if (unsigned(NSAM) < sample->size()/5)
	  {
	    sample->erase(sample->begin()+NSAM+1,sample->end());
	  }
      }
    if (max_chromosomes != NULL && sample->size() > *max_chromosomes)
      *max_chromosomes  += max_chromosomes_inc;

    //As we have scaled time in units of 4Nr generations, we pass theta
    //to the mutation function.  If we had used the following code to
    //generate times:
    //	double tcoal = expo(2./rcoal);
    //	double trec = expo(2./rrec);
    //Then time would be scaled in units of 2Ne generations, 
    //and if our simulation considers theta=4Neu, we'd pass theta/2
    //to the infinite_sites routine

    SimData gametes_obj = infinite_sites_sim_data(poiss,uni,
						  nsites,*sample_history,theta);
    return gametes_obj;
  }
}
