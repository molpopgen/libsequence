// -*- C++ -*-
/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

*/


/*! 
  \file Correlations.tcc
  @short implementation of functions declared in Correlations.hpp
*/
#ifndef __CORRELATIONS_TCC__
#define __CORRELATIONS_TCC__

#include <Sequence/Portability/random_shuffle.hpp>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <limits>

namespace Sequence 
{
  /**
     @ingroup stats
     Calculates Pearson's product-moment correlation between two containers, \c x and \c y.
     \f[r_{1,2}=\frac{y_1\times y_2}{\sqrt{\sum_{i=1}^{i=n} y_1^2 \times \sum \sum_{i=1}^{i=n} y_2^2}}\f]
     From Sokal and Rohlf. Biometry, 3rd ed. 1995, Freeman Press, page  571.
     @param beg_x pointer to beginning of \c x
     @param end_x pointer to one past the end of \c x
     @param beg_y pointer to beginning of \c y
     @note the square of the return value is the well-know coefficient of determination, r^2
  */
  template<typename iter1, typename iter2>
  typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
			  typename std::iterator_traits<iter2>::value_type>::type
  ProductMoment::operator()(iter1 beg_x, iter1 end_x, iter2 beg_y) const
  {
    typedef typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
      typename std::iterator_traits<iter2>::value_type>::type rtype;

    if (beg_x >= end_x) return std::numeric_limits<rtype>::min();

    unsigned nsam = 0;
    rtype x_bar=rtype();
    rtype y_bar=rtype();
    rtype x_squared=rtype();
    rtype y_squared=rtype();
    rtype xy=rtype();
    
    for ( ; beg_x != end_x  ; ++beg_x,++beg_y)
      {
	x_bar += *beg_x;
	x_squared += (*beg_x)*(*beg_x);
	y_bar += *beg_y;
	y_squared += (*beg_y)*(*beg_y);	
	xy += (*beg_x)*(*beg_y);
	++nsam;
      }
    x_bar /= rtype(nsam);
    y_bar /= rtype(nsam);
    x_squared /= rtype(nsam);
    y_squared /= rtype(nsam);
    xy -= x_bar*rtype(nsam)*y_bar;
    xy -= y_bar*rtype(nsam)*x_bar;
    xy += rtype(nsam)*x_bar*y_bar;
    rtype sum_x_sq = rtype(nsam)*(x_squared-x_bar*x_bar);
    rtype sum_y_sq = rtype(nsam)*(y_squared-y_bar*y_bar);
    return(xy/std::pow((sum_x_sq*sum_y_sq),0.5));
  }

  template<typename iter1, typename iter2>
  typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
			  typename std::iterator_traits<iter2>::value_type>::type
  SpearmansRank::operator()(iter1 beg_x, iter1 end_x, iter2 beg_y) const
  {
    typedef typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
      typename std::iterator_traits<iter2>::value_type>::type rtype;
  
  if(beg_x>=end_x) return std::numeric_limits<rtype>::min();
  
  typedef typename std::iterator_traits<iter1>::value_type t1;
  typedef typename std::iterator_traits<iter2>::value_type t2;
  typedef typename std::iterator_traits<iter1>::difference_type dtype1;
  typedef typename std::iterator_traits<iter2>::difference_type dtype2;
  dtype1 d = end_x-beg_x;
  std::vector<t1> xranks;
  std::vector<t2> yranks;
  //get the unique elements in each range
  std::unique_copy(beg_x,end_x,std::back_inserter(xranks));
  std::unique_copy(beg_y,beg_y+size_t(d),std::back_inserter(yranks));
  //sort them
  std::sort(xranks.begin(),xranks.end());
  std::sort(yranks.begin(),yranks.end());

  //Recode the data using a map<value,rank>,
  //as it's a quick way to do the rank lookups
  std::map<t1,unsigned> mx;
  std::map<t2,unsigned> my;
  unsigned rank=0;
  for (typename std::vector<t1>::const_iterator i =
	 xranks.begin() ;
       i != xranks.end() ;
       ++i)
    {
      mx[*i] = rank++;
    }
  rank=0;
  for (typename std::vector<t2>::const_iterator i =
	 yranks.begin() ;
       i != yranks.end() ;
       ++i)
    {
      my[*i] = rank++;
    }
  //calculate Pearson's rho on the ranks
  unsigned nsam = 0;
  rtype x_bar=rtype();
  rtype y_bar=rtype();
  rtype x_squared=rtype();
  rtype y_squared=rtype();
  rtype xy=rtype();
  typename std::vector<t1>::const_iterator a;
  typename std::vector<t2>::const_iterator b;
  rtype ra,rb;
  for( ; beg_x != end_x ; ++beg_x,++beg_y )
    {
      //find where each element is in the ranks vector
      ra = mx[*beg_x];
      rb = my[*beg_y];
      x_bar += ra;
      x_squared += ra*ra;
      y_bar += rb;
      y_squared += rb*rb;
      xy += ra*rb;
      ++nsam;
    }
    x_bar /= rtype(nsam);
    y_bar /= rtype(nsam);
    x_squared /= rtype(nsam);
    y_squared /= rtype(nsam);
    xy -= x_bar*rtype(nsam)*y_bar;
    xy -= y_bar*rtype(nsam)*x_bar;
    xy += rtype(nsam)*x_bar*y_bar;
    rtype sum_x_sq = rtype(nsam)*(x_squared-x_bar*x_bar);
    rtype sum_y_sq = rtype(nsam)*(y_squared-y_bar*y_bar);
    return(xy/std::pow((sum_x_sq*sum_y_sq),0.5));
  }
   
  template<typename iter1, typename iter2,  
	   typename correlation_type,
	   typename comparison_function,
	   typename UniformIntGenerator>
  typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
			  typename std::iterator_traits<iter2>::value_type>::type
  PermuteCorrelation_details(iter1 beg_x, iter1 end_x, iter2 beg_y,
			     const correlation_type &  corr,
			     const comparison_function & comp,
			     UniformIntGenerator & rand,
			     const unsigned & NPERM)
  {
    typedef typename  ensureFloating<typename std::iterator_traits<iter1>::value_type,
      typename std::iterator_traits<iter2>::value_type>::type rtype;
    typedef typename std::iterator_traits<iter1>::value_type type1;

    rtype _obs = corr(beg_x,end_x,beg_y);
    unsigned _prob=0;

    type1 *copy_x = new type1[end_x-beg_x+1];
    std::copy(beg_x,end_x,copy_x);

    for(unsigned i = 0 ; i < NPERM ; ++i)
      {
	Sequence::random_shuffle(copy_x,copy_x+(end_x-beg_x),rand);
	if ( comp(std::fabs(corr(copy_x,copy_x+(end_x-beg_x),beg_y)),
		  std::fabs(_obs)) == true )
	  {
	    ++_prob;
	  }
      }
    delete [] copy_x;
    return rtype(_prob)/rtype(NPERM);
  }

  template<typename iter1, typename iter2,  
	   typename correlation_type,
	   typename comparison_function,
	   typename UniformIntGenerator>
  typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
			  typename std::iterator_traits<iter2>::value_type>::type
  PermuteCorrelation(iter1 beg_x, iter1 end_x, iter2 beg_y,
		     const correlation_type &  corr,
		     const comparison_function & comp,
		     UniformIntGenerator & rand,
		     const unsigned & NPERM)
  /*! 
    Obtain the p-value of a correlation coefficient by permutation.  This function
    can be used to get 1- or 2- tailed p-values by using different comparison_function
    objects.  For example, using std::greater_equal<double> will returned the 1-tailed 
    probability of observing a correlation >= the observed value.
    @param beg_x pointer to the beginning of the range of the 1st vector
    @param end_x pointer to the end of the range of the 1st vector
    @param beg_y pointer to the beginning of the range of the 2nd vector
    @param corr a function object to calculate the correlation statistic (i.e. ProductMoment)
    @param comp a comparison function
    @param NPERM number of permutations to do
    @param rand a function returning a random integer (must be compatible with std::random_shuffle)
    @note This function keeps the order of the 2 containers intact.
  */
  {
    return PermuteCorrelation_details(beg_x,end_x,beg_y,corr,comp,rand,NPERM);
  }

  template<typename iter1, typename iter2,  
	   typename correlation_type,
	   typename comparison_function,
	   typename UniformIntGenerator>
  typename ensureFloating<typename std::iterator_traits<iter1>::value_type,
			  typename std::iterator_traits<iter2>::value_type>::type
  PermuteCorrelation(iter1 beg_x, iter1 end_x, iter2 beg_y,
		     const correlation_type &  corr,
		     const comparison_function & comp,
		     const UniformIntGenerator & rand,
		     const unsigned & NPERM)
  /*! 
    Obtain the p-value of a correlation coefficient by permutation.  This function
    can be used to get 1- or 2- tailed p-values by using different comparison_function
    objects.  For example, using std::greater_equal<double> will returned the 1-tailed 
    probability of observing a correlation >= the observed value.
    @param beg_x pointer to the beginning of the range of the 1st vector
    @param end_x pointer to the end of the range of the 1st vector
    @param beg_y pointer to the beginning of the range of the 2nd vector
    @param corr a function object to calculate the correlation statistic (i.e. ProductMoment)
    @param comp a comparison function
    @param NPERM number of permutations to do
    @param rand a function returning a random integer (must be compatible with std::random_shuffle)
    @note This function keeps the order of the 2 containers intact.
  */
  {
    return PermuteCorrelation_details(beg_x,end_x,beg_y,corr,comp,rand,NPERM);
  }

}//namespace Sequence
#endif
