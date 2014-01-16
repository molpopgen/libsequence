#ifndef __SEQUENCE_RNG_GSL_RNG_WRAPPERS_HPP__
#define __SEQUENCE_RNG_GSL_RNG_WRAPPERS_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace Sequence
{
  class gsl_uniform01
  {
  private:
    gsl_rng * __r;
  public:
    gsl_uniform01( gsl_rng * r ) : __r(r)
    {
    }
    inline double operator()(void) const
    {
      return gsl_rng_uniform(__r);
    }
  };

  class gsl_uniform
  {
  private:
    gsl_rng * __r;
  public:
    gsl_uniform( gsl_rng * r ) : __r(r)
    {
    }
    inline double operator()(double a, double b) const
    {
      return gsl_ran_flat(__r,a,b);
    }
  };

  class gsl_exponential
  {
  private:
    gsl_rng * __r;
  public:
    gsl_exponential( gsl_rng * r ) : __r(r)
    {
    }
    inline double operator()(const double & mean) const
    {
      return gsl_ran_exponential(__r,mean);
    }
  };

  class gsl_poisson
  {
  private:
    gsl_rng * __r;
  public:
    gsl_poisson( gsl_rng * r ) : __r(r)
    {
    }
    inline int operator()(const double & mean) const
    {
      return gsl_ran_poisson(__r,mean);
    }
  };
}
#endif

