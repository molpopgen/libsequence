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

// -*- C++ -*-
#ifndef __DESCRIPTIVE_STATS_TCC__
#define __DESCRIPTIVE_STATS_TCC__
#include <numeric>
#include <functional>
#include <iterator>

/*! \file descriptiveStats.tcc
  \short implementations for file descriptiveStats.hpp
*/

namespace Sequence
{
  template<typename iterator> 
  double mean(iterator beg, iterator end)
  {
    double _result(0);
    double _nsam(end-beg);
    for( ; beg != end ; beg++)
      {
	_result += *beg;
      }
    return double(_result)/double(_nsam);
  }

  template<typename iterator> 
  double variance(iterator beg, iterator end)
  {
    double _mean(0);
    double _xsq(0);
    double _nsam(end-beg);
    for( ; beg != end ; beg++)
      {
	_mean += *beg;
	_xsq += (*beg)*(*beg);
      }
    return double(_xsq)/double(_nsam-1) -
      (double(_mean)/double(_nsam))*(double(_mean)/double(_nsam-1));
  }

  template<typename ForwardIterator>
  std::pair<double,double>
  meanAndVar(ForwardIterator beg,
	     ForwardIterator end)
  {
    Sums<double> __s = std::accumulate(beg,end,Sums<double>());
    return std::make_pair(__s.mean(),__s.variance());
  }

  template<typename T>
  Sums<T>::Sums() : __sum( T() ), __sumsq( T() ), __n(0u)
  {
  }

  template<typename T>
  Sums<T> & Sums<T>::operator+=(const T & t)
  {
    __sum += t;
    __sumsq += t*t;
    __n++;
    return *this;
  }

  template<typename T>
  Sums<T> & Sums<T>::operator+=(const Sums<T> & s)
  {
    __sum += s.__sum;
    __sumsq += s.__sumsq;
    __n += s.__n;
    return *this;
  }

  template<typename T>
  const Sums<T> operator+(const Sums<T> & lhs,const Sums<T> & rhs)
  {
    return Sums<T>(lhs) += rhs;
  }

  template<typename T>
  const Sums<T> operator+(const Sums<T> & lhs,const T & rhs)
  {
    return Sums<T>(lhs) += rhs;
  }
  
  template<typename T>
  const T & Sums<T>::sum() const
  {
    return __sum;
  }

  template<typename T>
  const T & Sums<T>::sumSquares() const
  {
    return __sumsq;
  }

  template<typename T>
  typename Sums<T>::floating_type Sums<T>::mean() const
  {
    return floating_type(__sum)/floating_type(__n);
  }

  template<typename T>
  typename Sums<T>::floating_type Sums<T>::variance() const
  {
    return ( floating_type(__sumsq)/floating_type(__n-1) -
	     floating_type(__sum*__sum)/floating_type(__n*(__n-1)) );
  }
}
#endif
