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

#ifndef SIMDATA_H
#define SIMDATA_H
/*! \file SimData.hpp
  @brief Declaration of Sequence::SimData, a class representing polymorphism data
  from coalescent simulations generated under an infinite-sites scheme
*/

/*! \class Sequence::SimData Sequence/SimData.hpp
  \ingroup polytables
  \ingroup coalescent
 A class for reading in data from ms. ms is Dick
 Hudson's coalescent simulation program, available from
 http://home.uchicago.edu/~rhudson1/.  This class reads in
 the data one record at a time.  As of the time this document
 was written, the output of Hudon's program for one run looks like:\n
 \n
 //\n
 segsites: 10\n
 positions: 0.2512 0.3212 0.3449 0.4386 0.4571 0.4842 0.5745 0.6333 0.7042 0.9928 \n
 0100010010\n
 0000000100\n
 0100010011\n
 0000000000\n
 0100010010\n
 0101010010\n
 0100010011\n
 0110010010\n
 1100111010\n
 0000000000\n
 \n
 Look at the example program tajd.cc for an example of how to pipe
 simulated data into a C++ program using these libraries.
\note Be sure that you notice that the positions of the segregating sites
fall in the interval (0,1].
 @author Kevin Thornton
 @short Data from coalescent simulations
*/

/*! \example msstats.cc
 */

#include <Sequence/PolyTable.hpp>
#include <cstdio>
namespace Sequence
  {
  class SimData:public PolyTable
    {
    private:
      size_t totsam;
    public:
      explicit SimData (const size_t & nsam=0, const size_t & nsnps = 0);
      explicit SimData(double *pos, char **sample, int nsam, int S);
      explicit SimData(const std::vector<double> & pos, const std::vector<std::string> & data)
      {
        assign(&pos[0],unsigned(pos.size()),&data[0],unsigned(data.size()));
      }
      explicit SimData(const SimData::const_site_iterator beg, 
		       const SimData::const_site_iterator end) : PolyTable(beg,end)
      {
      }
      ~ SimData (void){}
      
      virtual std::istream & read (std::istream & s) ;
      virtual std::ostream & print(std::ostream &o) const;

      void Binary (bool haveOutgroup = false,
                   unsigned outgroup = 0,
                   bool strictInfSites = true)
      {
        //no need to do anything...
        return;
      }
      virtual int fromfile( FILE * openedfile );
      inline unsigned segsites (void) const
      /*!
        Returns the number of segregating sites in the 
        data block.
      */
        {
          return (this->numsites());
        }
    };

}
#endif
