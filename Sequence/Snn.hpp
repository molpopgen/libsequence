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

#ifndef __SEQUENCE_SNN_HPP__
#define __SEQUENCE_SNN_HPP__

#include <Sequence/Comparisons.hpp>
#include <Sequence/PolyTable.hpp>
#include <vector>
#include <utility>
#include <cassert>

namespace Sequence
{
  double Snn_statistic( const unsigned individuals[],
			const std::vector< std::vector<double> > & dkj,
			const unsigned config[],
			const size_t & npop,
			const unsigned & nsam );

  template< typename uniform_int_generator >
  std::pair<double,double>
  Snn_test(const PolyTable & snpTable,
	   const unsigned config[],
	   const size_t & npop,
	   uniform_int_generator & uni_int,
	   const unsigned & nperms = 10000);

  template< typename uniform_int_generator >
  std::vector< std::vector<double> >
  Snn_test_pairwise(const PolyTable & snpTable,
		    const unsigned config[],
		    const size_t & npop,
		    uniform_int_generator & uni_int,
		    const unsigned & nperms = 10000);

  template< typename uniform_int_generator >
  std::pair<double,double>
  Snn_test(const PolyTable & snpTable,
	   const unsigned config[],
	   const size_t & npop,
	   const uniform_int_generator & uni_int,
	   const unsigned & nperms = 10000);

  template< typename uniform_int_generator >
  std::vector< std::vector<double> >
  Snn_test_pairwise(const PolyTable & snpTable,
		    const unsigned config[],
		    const size_t & npop,
		    const uniform_int_generator & uni_int,
		    const unsigned & nperms = 10000);
}
#endif
#include <Sequence/bits/Snn.tcc>
