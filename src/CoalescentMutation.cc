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

#include <Sequence/Coalescent/Mutation.hpp>

namespace Sequence
{
  void output_gametes(FILE * fp,
		      const unsigned & segsites,
		      const unsigned & nsam,
		      const gamete_storage_type & gametes)
  /*!
    @brief Write an object of type gamete_storage type to a C-style file stream
    This function is used when you need to output simulated gametes using a
    method faster than the operator<< for class SimData.
    \param fp pointer to an open C-style output stream
    \param segsites the number of segregating sites in \a gametes
    \param nsam the number of individuals in \a gametes
    \param gametes the simulated sample.  Must be allocated to hold at least 
    \a segsites positions, and \a nsam strings of length \a segsites
  */
  {
    fprintf(fp,"//\n");
    if ( segsites > 0 )
      {
	fprintf(fp,"segsites: %d\npositions: ",segsites);
	for(unsigned i=0;i<segsites;++i)
	  {
	    fprintf(fp,"%f ",gametes.first[i]);
	  }
	fprintf(fp,"\n");
	for(unsigned i=0;i<nsam;++i)
	  {
	    for(unsigned j=0;j<segsites;++j)
	      {
		fprintf(fp,"%c",gametes.second[i][j]);
	      }
	    fprintf(fp,"\n");
	  }
      }
    else
      {
	fprintf(fp,"segsites: 0\n");
      }
  }
}
