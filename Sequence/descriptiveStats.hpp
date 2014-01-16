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

#ifndef __DESCRIPTIVE_STATS_H__
#define __DESCRIPTIVE_STATS_H__
#include <Sequence/ensureFloating.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <utility>

/*! \file descriptiveStats.hpp 
  \short classes and functions for calculating descriptive statistics of lists of random variables
*/


/*! \defgroup stats Statistics
  @short routines for statistics
*/
namespace Sequence
{
  /*!
    \param beg an iterator
    \param end an iterator
    \return the mean of the range
    \ingroup stats
  */
  template<typename iterator> 
  double mean(iterator beg, iterator end);

  /*!
    \param beg an iterator
    \param end an iterator
    \return the variance of the range
    \ingroup stats
  */
  template<typename iterator> 
  double variance(iterator beg, iterator end);

  /*!
    A function to calculate the mean and variance of the values
    stored in a container.  The rationale is that when both
    the mean and the variance (an sum of squares) are needed,
    it is more efficient to calculate them together, because
    you only go over the data once.
  */
  template<typename ForwardIterator> 
  std::pair<double,double> meanAndVar(ForwardIterator beg,
				      ForwardIterator end);

  /*!
    \class Sums Sequence/descriptiveStats.hpp
    A class to keep track of the sum, and sum of squares,
    of values.  I wrote it to make it easier to calculate
    the means and variances.

    Examples:
    \code
    //example illustrates use of += and using std::accumulate
    Sequence::Sums<unsigned> s1;
    std::vector<unsigned> vd;
    int i = 0;
    while( i++ < 100 )
    {
    unsigned r = random;
    s1 += r;
    vd.push_back(r);
    }
    Sequence::Sums<double> s2 = std::accumulate(vd.begin(),vd.end(),Sequence::Sums<double>());
    std::cout << s1.mean() << '\t' << s2.mean() << '\n'
    << s1.variance() << '\t' << s2.variance() << '\n';
    \endcode
    \note Can be used with std::accumulate
  */
  template< typename T >
  class Sums
  {
  public:
    typedef typename ensureFloating<T,T>::type floating_type;
    //T must be convertible to floating_type via a typecast
    BOOST_STATIC_ASSERT( (boost::is_convertible<T,floating_type>::value) );
  private:
    T __sum,__sumsq;
    unsigned __n;
  public:
    Sums();
    Sums<T> & operator+=(const T &);
    Sums<T> & operator+=(const Sums<T> &);
    const T & sum() const;
    const T & sumSquares() const;
    floating_type mean() const;
    floating_type variance() const;
  };

  template<typename T>
  const Sums<T> operator+(const Sums<T> & lhs,const Sums<T> & rhs);

  template<typename T>
  const Sums<T> operator+(const Sums<T> & lhs,const T & rhs);

}//namespace Sequence
#include <Sequence/bits/descriptiveStats.tcc>
#endif
