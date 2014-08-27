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

#ifndef __CORRELATIONS_H__
#define __CORRELATIONS_H__

/*! \file Correlations.hpp
  
  Generic algorithms for calculating and testing the significance
  of correlation coefficients
*/
/*! \example correlations.cc */
#include <iterator>
#include <functional>

namespace Sequence
{
  struct ProductMoment
  /*!
    Function object to calculate Pearson's product-moment correlation
    coefficient
    @ingroup stats
    @short Pearson's product-moment correlation
  */
  {
    template<typename iter1, typename iter2>
    double operator()(iter1 beg_x,iter1 end_x, iter2 beg_y) const;
  };

  struct SpearmansRank
  /*!
    Function object to calculate Spearman's rank correlation
    @ingroup stats
    @short Spearman's rank correlation
  */
  {
    template<typename iter1,typename iter2>
    double operator()(iter1 beg_x,iter1 end_x, iter2 beg_y) const;
  };

  template<typename iter1, typename iter2,  
	   typename correlation_type,
	   typename comparison_function,
	   typename shuffler>
  double PermuteCorrelation(iter1 beg_x, iter1 end_x, iter2 beg_y,
			    const correlation_type &  c,
			    const comparison_function & comp,
			    shuffler & s,
			    const unsigned & NPERM=10000);

  
  template<typename iter1, typename iter2,  
	   typename correlation_type,
	   typename comparison_function,
	   typename shuffler>
  double PermuteCorrelation(iter1 beg_x, iter1 end_x, iter2 beg_y,
			    const correlation_type &  c,
			    const comparison_function & comp,
			    const shuffler & s,
			    const unsigned & NPERM=10000);
  
}//namespace Sequence
#include <Sequence/bits/Correlations.tcc>
#endif
