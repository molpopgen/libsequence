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

#include <Sequence/Snn.hpp>
namespace Sequence
{
  double Snn_statistic( const unsigned individuals[],
			const std::vector< std::vector<double> > & dkj,
			const unsigned config[],
			const size_t & npop,
			const unsigned & nsam)
{
  /*
    notation for variables follows Hudson's paper
  */
  double snn = 0.;

  //store the d_kj for the whole sample
  double * d_kj = new double[nsam-1];

  //store d_kj for within-population comparisons
  std::vector<double> d_kj_win;
  for(unsigned k=0; k<nsam ; ++k)
    {
      //find out which pop ind k is in
      unsigned pop = 0,ttl=0;
      while (pop < npop)
	{
	  ttl += config[pop];
	  if (k < ttl)
	    break;
	  pop++;
	}
      d_kj_win.clear();
      for(unsigned j = 0,dummy=0; j < nsam ; ++j)
	{
	  if (k!=j)
	    {
	      unsigned a=*(individuals+j),b=*(individuals+k);
	      if (a>b)
		std::swap(a,b);

	      double ndiffs = dkj[a][b];
	      d_kj[dummy++] = ndiffs;
	      //figure out what pop j is in;
	      unsigned pop_j=0,ttl=0;
	      while (pop_j < npop)
		{
		  ttl += config[pop_j];
		  if (j < ttl)
		    break;
		  pop_j++;
		}
	      if (pop==pop_j)
		d_kj_win.push_back(ndiffs);
	    }
	}
      //Calculate T_k
      double min = d_kj[0];
      for (unsigned j = 1 ; j < nsam-1 ; ++j)
	if (d_kj[j] < min) min = d_kj[j];
      
      std::ptrdiff_t T_k = std::count(d_kj,d_kj+(nsam-1),min);

      //Calculate M_k
      std::ptrdiff_t M_k = std::count(d_kj_win.begin(),
				d_kj_win.end(),min);
      snn += double(M_k)/double(T_k);
    }
  delete [] d_kj;
  return snn/double(nsam);
}

}
