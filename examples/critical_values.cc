#include <Sequence/Crit.hpp>
#include <Sequence/descriptiveStats.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>

using namespace Sequence;
using namespace std;

int main()
{
  //fill a vector<double> with 1000 random values
  std::mt19937 generator(std::time(0));
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  vector<double> x;
  for(unsigned i = 0 ; i < 100000 ; ++i)
    {
      x.push_back( distribution(generator) );
    }

  //sort the list
  sort(x.begin(),x.end());

  //get upper 95% critical value
  std::pair< std::vector<double>::iterator, double > p = 
    upperCrit()(x.begin(),x.end());
  cout << *(p.first) << '\t' << p.second << endl;

  //get lower 5% critical value
  p = lowerCrit()(x.begin(),x.end());
  cout << *(p.first) << '\t' << p.second << endl;

  //examples of descriptive stats
  std::pair<double,double> desc = meanAndVar(x.begin(),x.end());
  cout << desc.first << '\t' << desc.second << endl;
  cout << mean(x.begin(),x.end()) << '\t'
       << variance(x.begin(),x.end()) << endl;
}
