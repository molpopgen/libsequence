/*! \include  bottleneck.cc
*/

/*
  Simulates 10000 samples from a bottleneck model where
  the population recovers from the bottleneck according to an 
  exponential growth model.

  The program does not output the samples, but instead
  calculates the number of segregating sites, pi,
  Tajima's D, and haplotype diversity for each sample.

  With respect to Dick Hudson's program "ms", the equivalent parameters for ms
  would be:

  ms 10 10000 -t 10 -r 10 1000 -eG 0.1 8.047 -eG 0.3 0 -eN 0.3 1
*/

#include <Sequence/Coalescent/DemographicModels.hpp>
#include <Sequence/Coalescent/Initialize.hpp>
#include <Sequence/Coalescent/Mutation.hpp>
#include <Sequence/PolySIM.hpp>
#include <iostream>
#include <boost/bind.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>

int main(int argc, char **argv)
{
  gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r,std::time(0));

  const unsigned n = 10;
  std::vector<Sequence::chromosome> sample = Sequence::init_sample(std::vector<int>(1,n),1000);
  Sequence::marginal imarg = Sequence::init_marginal(n);
  for(unsigned i=0;i<10000;++i)
    {
      //simulate the ancestral recombination graph for the sample
      Sequence::arg hist = Sequence::bottleneck(boost::bind(gsl_ran_flat,r,_1,_2),
						boost::bind(gsl_rng_uniform,r),
						boost::bind(gsl_ran_exponential,r,_1),
						sample,imarg,0.1,0.2,0.2,10.,true,1.);
      //Apply mutations according to the infinitely-many sites scheme
      Sequence::SimData d = Sequence::infinite_sites_sim_data(boost::bind(gsl_ran_poisson,r,_1),
							      boost::bind(gsl_ran_flat,r,_1,_2),
							      1000,hist,10.);
      Sequence::PolySIM ad(&d);
      std::cout << d.numsites() << ' ' << ad.ThetaPi() << ' '
		<< ad.TajimasD() << ' ' << ad.DandVH() << '\n';
    }
}
