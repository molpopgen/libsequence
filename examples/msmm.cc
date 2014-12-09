/*! \include ms--.cc
*/

#include <Sequence/Coalescent/NeutralSample.hpp>
#include <Sequence/Coalescent/Initialize.hpp>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <random>

/*
  These integers are declared extern in <Sequence/Coalescent/Mutation.hpp>,
  and are used to make storage of gametes more efficient during the simulation.
  You need to set their value here.  In practice, the values below work well.
*/
namespace Sequence{
  MAX_SEG_T MAX_SEGSITES = 200;
  MAX_SEG_T MAX_SEGS_INC = 100;
}

int main(int argc, char **argv)
{
  //get basic parameters from command line
  const int nsam = atoi(argv[1]);
  int howmany = atoi(argv[2]);
  const double theta = atof(argv[3]);
  const double rho = atof(argv[4]);
  const int nsites = atoi(argv[5]);

  //Set up our random number distribution functions
  unsigned seed = std::time(0);
  std::mt19937 generator(seed);
  //Make our RNG types via lambda expressions
  auto poiss = [&generator](const double & mean){ return std::poisson_distribution<int>(mean)(generator); };
  auto expo = [&generator](const double & mean){ return std::exponential_distribution<double>(1./mean)(generator); };
  auto uni01 = [&generator](){ return std::uniform_real_distribution<double>(0.,1.)(generator); };
  auto uni = [&generator](const double & a, const double & b){ return std::uniform_real_distribution<double>(a,b)(generator); };


  //initialize a vector of chromosomes
  //There will be 1 population containing nsam chromosomes with nsites each
  //sites are labelled 0 to nsites-1
  std::vector<Sequence::chromosome> initialized_sample = Sequence::init_sample( std::vector<int>(1,nsam),nsites );
  //initialize a marginal tree
  Sequence::marginal initialized_marginal = Sequence::init_marginal(nsam);

  //These value can be used by Sequence::neutral_sample
  //to fine-tune memory allocation in the simulation
  unsigned MAXCH=0;
  const unsigned MAXCH_INC=50;
  double sum=0.;
  
  for(int i=1;i<nsam;++i)
    sum+= 1./double(i);
  MAXCH=5*unsigned(rho*sum);

  for(int i=0;i<argc;++i)
    {
      std::cout << argv[i] << ' ';
    }
  std::cout << '\n';

  while(howmany--)
    {
      //copy construct to avoid having to re-init each time
      std::vector<Sequence::chromosome> sample;
      sample.reserve(MAXCH); //use MAXCH to reserve a good amount of contiguous memory
      std::copy(initialized_sample.begin(),
		initialized_sample.end(),
		std::back_inserter(sample));
      Sequence::arg sample_history(1,initialized_marginal);

      //simulate a sample from the standard neutral model
      //of a large Wright-Fisher population with infinite-sites 
      //mutation
      Sequence::SimData gametes = Sequence::neutral_sample(uni,uni01,expo,poiss,
							   theta,rho,nsites,nsam,
							   &sample,
							   &sample_history,
							   &MAXCH,
							   MAXCH_INC);
      //print it
      std::cout << gametes << '\n';
    }
}
