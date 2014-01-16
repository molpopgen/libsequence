/*! \include msbeta.cc
*/

#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/RNG/gsl_rng_wrappers.hpp>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <gsl/gsl_cdf.h>

/*
  These integers are declared extern in <Sequence/Coalescent/Mutation.hpp>,
  and are used to make storage of gametes more efficient during the simulation.
  You need to set their value here.  In practice, the values below work well.
*/
namespace Sequence{
  int MAX_SEGSITES = 200;
  int MAX_SEGS_INC = 100;
}

int main(int argc, char **argv)
{
  //get basic parameters from command line
  const int nsam = atoi(argv[1]);
  int howmany = atoi(argv[2]);
  const double theta = atof(argv[3]);
  const double rho = atof(argv[4]);
  const int nsites = atoi(argv[5]);
  const double a = atof(argv[6]);
  const double b = atof(argv[7]);

  //initialize a mersenne twister and seed it with system time
  //These types are defined in the GNU Scientific Library (GSL)
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,time(0));

  //Wrapper types for gsl functions declared in
  //<Sequence/RNG/gsl_rng_wrappers.hpp>
  //These all declare operator()
  Sequence::gsl_uniform01 uni01(r); //takes no arguments
  Sequence::gsl_uniform uni(r);     //takes two doubles a and b, and returns a U[a,b)
  Sequence::gsl_exponential expo(r);//takes the mean as an argument
  Sequence::gsl_poisson poiss(r);   //takes the mean as an argument

  //initialize a vector of chromosomes
  //There will be 1 population containing nsam chromosomes with nsites each
  //sites are labelled 0 to nsites-1
  std::vector<Sequence::chromosome> initialized_sample = Sequence::init_sample( std::vector<int>(1,nsam),nsites );
  //initialize a marginal tree
  Sequence::marginal initialized_marginal = Sequence::init_marginal(nsam);

  //defined in Sequence/Coalescent/Mutation.hpp
  Sequence::gamete_storage_type gametes( std::vector<double>(Sequence::MAX_SEGSITES,0.),
					 std::vector<std::string>(nsam,
								  std::string(Sequence::MAX_SEGSITES,'0')) );
  int segsites;

  //These value can be used by Sequence::neutral_sample
  //to fine-tune memory allocation in the simulation
  unsigned MAXCH=0;
  double sum=0.;
  
  for(int i=1;i<nsam;++i)
    sum+= 1./double(i);
  MAXCH=5*unsigned(rho*sum);

  for(int i=0;i<argc;++i)
    {
      std::cout << argv[i] << ' ';
    }
  std::cout << '\n';

  //first, we need to make the pdf for our genetic map,
  //which is determined by a Beta(a,b) density * the nsites-1
  std::vector<double> discretized_beta(nsites-1,0.);
  for(int i = 0 ; i < nsites-1 ; ++i)
    {
      discretized_beta[i] = gsl_cdf_beta_P( double(i+1)/double(nsites-1),a,b ) -
	gsl_cdf_beta_P( double(i)/double(nsites-1),a,b );
    }

  while(howmany--)
    {
      std::vector<Sequence::chromosome> sample;
      sample.reserve(MAXCH);
      std::copy(initialized_sample.begin(),
		initialized_sample.end(),
		std::back_inserter(sample));
      int NSAM = nsam;
      Sequence::arg sample_history(1,initialized_marginal);
      int nlinks = nsam*(nsites-1);
      double t = 0.;

      //generate samples under a WF model, scaling time in units of 4Ne generations
      while(NSAM>1)
	{
	      double rcoal = double(NSAM*(NSAM-1));
	      double rrec = 0.;
	      double reclen = 0.;
	      std::vector<double> precs;
	      
	      //since the recombination map is non-uniform,
	      //the contribution of each (segment of a ) chromosome
	      //is not equal, so we need to sum up rrec
	      //by integrating the beta density over the length of 
	      //each chromosome
	      if(rho>0)
		{
		  for(int i=0;i<NSAM;++i)
		    {
		      int beg = sample[i].first();
		      int end = sample[i].last();
		      double cumulative = std::accumulate(discretized_beta.begin()+beg,
							  discretized_beta.begin()+end,0.);
		      rrec += cumulative*(double(end-beg));
		      reclen += cumulative;
		      precs.push_back(cumulative);
		    }
		}

	      //scale time in units of 4Ne generations
	      double tcoal = expo(1./rcoal);
	      double trec = expo(1./rrec);
	      if ( trec < tcoal ) //crossover event
		{
		  t+=trec;
		  //pick the recombinant chromosome, and the crossover position,
		  //based on our genetic map
		  std::pair<int,int> pos_rec = Sequence::pick_spot(uni01,
								   reclen,
								   precs,
								   sample.begin(),
								   NSAM,&discretized_beta[0]);
		  
		  nlinks -= Sequence::crossover(NSAM,pos_rec.first,pos_rec.second,
						&sample,&sample_history);
		  NSAM++;
		}
	      else //common ancestor event
		{
		  t+=tcoal;
		  std::pair<int,int> two = Sequence::pick2(uni,NSAM);
		  NSAM -= Sequence::coalesce(t,nsam,NSAM,two.first,two.second,nsites,
					     &nlinks,&sample,&sample_history);
		}
	      if (NSAM < sample.size()/5)
		{
		  sample.erase(sample.begin()+NSAM+1,sample.end());
		}
	    }

      //Uncomment lines below to print out the beginnings of all marginal trees to stderr
      //if you redirect these into a file, and plot them (dividing by nsites-1), 
      //the distribution will be the Beta(a,b) provided that nsites is reasonably large.
      //Note that there will be howmany values equal to 0, as each sample must have
      //at least 1 marginal beginning at 0.
      //for(Sequence::arg::iterator ti = sample_history.begin();ti!=sample_history.end();++ti)
      //{
      //std::cerr << ti->beg << '\n';
      //}
      segsites = Sequence::infinite_sites(poiss,uni,&gametes,
					  nsites,sample_history,theta);
      Sequence::output_gametes(stdout,segsites,nsam,gametes);

      if(sample.size()>MAXCH) MAXCH=2*sample.size();
      sample.clear();
      sample_history.clear();
    }
}

