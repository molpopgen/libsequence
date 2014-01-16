#ifndef __RANDOM_SHUFFLE_FIXER_HPP__
#define __RANDOM_SHUFFLE_FIXER_HPP__

#include <functional>

/*! 
  \defgroup portability Portability
*/

namespace Sequence
{
  template<typename uni01>
  struct randomShuffleAdaptor : public std::binary_function<int,uni01,int>
				/*!
				  A functor to help adapt random number generators
				  to compile with std::random_shuffle.  This is also 
				  useful when random_shuffle won't compile when instantiated
				  without a random number generator (such as some cygwin versions).

				  The following shows example uses of this functor using two different
				  portable random number generation systems--BOOST (http://ww.boost.org),
				  and the gsl (http://sources.redhat.com/gsl).  The gsl example
				  is contrived, as one could make a more direct functor to adapt gsl 
				  routines, but you get the point.
				  \ingroup portability
				  \code
				  #include <Sequence/Portability/randomShuffleAdaptor.hpp>
				  #include <boost/random/mersenne_twister.hpp>
				  #include <boost/random/uniform_01.hpp>
				  #include <gsl/gsl_rng.h>
				  #include <gsl/gsl_randist.h>
				  #include <iterator>
				  #include <algorithm>
				  #include <vector>
				  #include <ctime>

				  typedef boost::mt19937 RNG;
				  typedef boost::uniform_01< RNG > _prob;

				  template<typename iterator>
				  void print(iterator beg, iterator end)
				  {
				  std::copy(beg,end,
				  std::ostream_iterator< typename std::iterator_traits<iterator>::value_type > (std::cout, " "));
				  std::cout << '\n';
				  }

				  struct randomU : public std::unary_function<void,double>
				  {
				  gsl_rng * __r;
				  explicit randomU(gsl_rng * r) : __r(r) {}
				  inline double operator()(  ) const
				  {
				  return gsl_ran_flat(__r,0.,1.);
				  }
				  };

				  int main(int argc, char **argv)
				  {
				  RNG generator;
				  _prob uni(generator);
				  std::vector<int> data;

				  gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
				  gsl_rng_set(r,time(0));

				  for(int i=0;i<10;++i)
				  {
				  data.push_back(i);
				  }

				  std::cout << "original data\n";
				  print(data.begin(),data.end());

				  randomShuffleAdaptor<_prob> ru(&uni);
				  std::random_shuffle(data.begin(),data.end(),ru);

				  std::cout << "shuffle with boost:\n";
				  print(data.begin(),data.end());

				  randomU rugsl(r);
				  randomShuffleAdaptor<randomU> ru2(&rugsl);
				  std::random_shuffle(data.begin(),data.end(),ru2);
				  std::cout << "shuffle with gsl:\n";
				  print(data.begin(),data.end());
				  }
				  \endcode
				*/
  {
    uni01 * uni;
    randomShuffleAdaptor( uni01 * u ) : uni(u) {}
    inline int operator()(const int & i) const
    {
      return int( (*uni)() * i);
    }
  };
}
#endif
