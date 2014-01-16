#include <Sequence/Correlations.hpp>
#include <Sequence/RNG/gsl_rng_wrappers.hpp>
#include <Sequence/Portability/randomShuffleAdaptor.hpp>
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>

int main()
{
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,std::time(0));

  std::vector<int> x(12);
  std::vector<double> y(12);
      
  x[0]=159;
  x[1]=179;
  x[2]=100;
  x[3]=45;
  x[4]=384;
  x[5]=230;
  x[6]=100;
  x[7]=320;
  x[8]=80;
  x[9]=220;
  x[10]=320;
  x[11]=210;
  y[0]=14.4;
  y[1]=15.2;
  y[2]=11.3;
  y[3]=2.5;
  y[4]=22.7;
  y[5]=14.9;
  y[6]=1.41;
  y[7]=15.81;
  y[8]=4.19;
  y[9]=15.39;
  y[10]=17.25;
  y[11]=9.52;

  std::cout << "original data:\n";
  std::copy(x.begin(),x.end(),std::ostream_iterator<int>(std::cout," "));
  std::cout << std::endl;
  std::copy(y.begin(),y.end(),std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
  std::cout << "Product-moment correlation: "
	    << Sequence::ProductMoment()( x.begin(),x.end(),y.begin() )
	    <<std::endl
	    << "Spearman's rank correlation: "
	    << Sequence::SpearmansRank()( x.begin(),x.end(),y.begin() )<<std::endl;
  Sequence::gsl_uniform01 uni01(r);
  Sequence::randomShuffleAdaptor<Sequence::gsl_uniform01> ru(&uni01);
  std::cout << "(all pvalues are one-tailed tests of observing a more extreme,\n"
	    << "i.e. larger correlation than the observe dvalue)\n"
	    << "pval of P-M correlation, seed with time: "
	    << Sequence::PermuteCorrelation(x.begin(),x.end(),y.begin(),
 					    Sequence::ProductMoment(),
					    std::greater_equal<double>(),ru)<<std::endl
	    << "pval of rank correlation, default seed: "
	    << Sequence::PermuteCorrelation(x.begin(),x.end(),
  					    y.begin(),
					    Sequence::SpearmansRank(),
					    std::greater_equal<double>(),ru) << std::endl
	    << "pval of rank correlation, seed with time: "
    	    << Sequence::PermuteCorrelation(x.begin(),x.end(),y.begin(),
 					    Sequence::SpearmansRank(),
					    std::greater_equal<double>(),ru)<<std::endl;
  std::cout << "original data are intact:\n";
  std::copy(x.begin(),x.end(),std::ostream_iterator<int>(std::cout," "));
  std::cout << std::endl;
  std::copy(y.begin(),y.end(),std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
}
