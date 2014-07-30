#include <Sequence/Correlations.hpp>
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <random>

//simple wrapper for us to use
void shuffle_it( std::vector<int>::iterator beg,
		 std::vector<int>::iterator end,
		 std::mt19937 & generator )
{
  std::shuffle(beg,end,generator);
}

int main()
{
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

  std::mt19937 generator(std::time(0));

  std::cout << "original data:\n";
  for( const auto _x : x )
    {
      std::cout << _x << ' ';
    }
  std::cout << std::endl;
  for ( const auto _y : y ) 
    {
      std::cerr << _y << ' ';
    }
  std::cout << std::endl;
  std::cout << "Product-moment correlation: "
	    << Sequence::ProductMoment()( x.begin(),x.end(),y.begin() )
	    <<std::endl
	    << "Spearman's rank correlation: "
	    << Sequence::SpearmansRank()( x.begin(),x.end(),y.begin() )<<std::endl;


  std::function<void(std::vector<int>::iterator,
		     std::vector<int>::iterator)> __x = [&generator](std::vector<int>::iterator  a,
								     std::vector<int>::iterator  b) { return shuffle_it(a,b,generator); };

  std::cout << "(all pvalues are one-tailed tests of observing a more extreme,\n"
	    << "i.e. larger correlation than the observe dvalue)\n"
	    << "pval of P-M correlation, seed with time: "
	    << Sequence::PermuteCorrelation(x.begin(),x.end(),y.begin(),
 					    Sequence::ProductMoment(),
					    std::greater_equal<double>(),
					    __x) << '\n'
	    << "pval of rank correlation, default seed: "
	    << Sequence::PermuteCorrelation(x.begin(),x.end(),
					    y.begin(),
					    Sequence::SpearmansRank(),
					    std::greater_equal<double>(),
					    __x) << std::endl
	    << "pval of rank correlation, seed with time: "
    	    << Sequence::PermuteCorrelation(x.begin(),x.end(),y.begin(),
 					    Sequence::SpearmansRank(),
					    std::greater_equal<double>(),
					    __x) <<std::endl;

  std::cout << "original data are intact:\n";
  for( const auto _x : x )
    {
      std::cout << _x << ' ';
    }
  std::cout << std::endl;
 for( const auto _y : y )
    {
      std::cout << _y << ' ';
    }
  std::cout << std::endl;
}
