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

// Code for the -*- C++ -*- namespace Sequence::PolySites template members
#include <Sequence/stateCounter.hpp>
#include <Sequence/SeqConstants.hpp>
#if (!__STRICT_ANSI__)
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Alignment.hpp>
#endif
#include <cctype>
namespace
{
  /*!
    It is often the case that alignment programs replace characters
    at a site with a period (.) if that character is the same as that
    found at that position in the first sequence in the alignment.
    This variable is used by Sequence::PolySites when turning aligned
    data into a polymorphism table.
  */
  const char IDENTICAL = '.';
}//anon namespace

namespace Sequence
{
   template<typename __DataType>
   PolySites::PolySites (const std::vector < __DataType >&alignment,
			 bool strictInfSites , bool ignoregaps,
			 bool skipMissing,  bool skipAdjSNP,
			 unsigned freqfilter):
    PolyTable()
    /*!
      This is the constructor if you are using "string-like" data, such as std::string, or
      Sequence::Fasta.  Note that the vector name is 
      aligment, and that means that every sequence had better be the same length!\n
      \n
      By default, there is no limit to how many characters can "segregate" at a
      variable position, although if there are more than 4, most biologists will
      start to worry.  There are, however, times when you may wish to onlu consider
      sites that have a total of 2 character states.  (NOTE: by two states, I mean
      including BOTH the ingroup and the outgroup sequence.) Setting strictInfSites
      to 1 will result in making a polymorphic sites object containing only sites with
      2 states.\n
      \param alignment vector of data
      \param strictInfSites if \c true, throw out all sites with > 2 states
      \param ignoregaps if \c true, do not count gapped sites as polymorphisms
      \param skipMissing if \c true, ignore ALL sites with missing data ('N')
      \param skipAdjSNP if \c does nothing.  a placeholder for a future feature
      \param freqfilter Defaults to 0. For a polymorphic site to be included in the
      final table, the minor allele count in the data (i.e. the number of times
      the minor allele occurs at that site) must be strictly greater than 
      \c freqfilter
      \note segsite positions are stored as positions (starting from 1)
      \warning when ignoregaps=false, this class does not do the right thing
    */
  {
    if (Alignment::IsAlignment(alignment) == true)
      {
        seqlen = alignment[0].length ();
        numseqs = alignment.size ();
        fillIt<__DataType> (alignment, strictInfSites, ignoregaps,skipMissing,freqfilter);
        if(skipAdjSNP == true)
	  {}
      }
  }

  template<typename __DataType>
  void
  PolySites::fillIt (const std::vector < __DataType >&alignment, bool strict, bool ignoregaps,
                     bool skipMissing,unsigned freqfilter)
  {
    std::vector<double> _positions;

    for(unsigned j =0 ; j < seqlen ; ++j)
      {
        stateCounter Counts;
        bool haveMissing = false;
        for(unsigned i = 0 ; i < numseqs ; ++i)
          {
            if (std::toupper(alignment[i][j]) == 'N'||
		(alignment[i][j]==IDENTICAL && std::toupper(alignment[0][j]) == 'N'))
              {
                haveMissing = true;
              }
            if (alignment[i][j] == IDENTICAL)
              {
                Counts(alignment[0][j]);
              }
            else
              {
                Counts(alignment[i][j]);
              }
          }
        if ( (ignoregaps == true && Counts.gap == 0) || ignoregaps == false )//skip gaps if necessary
          {
            if (skipMissing == false || (skipMissing == true && haveMissing == false))
              {
                unsigned numstates = Counts.nStates();
                unsigned minor_count = SEQMAXUNSIGNED;
                if ( freqfilter > 0 )
                  {
                    minor_count = (Counts.a > 0 && (Counts.a < minor_count || 
						    minor_count == SEQMAXUNSIGNED))
		      ? Counts.a : minor_count;
                    minor_count = (Counts.g > 0 && (Counts.g < minor_count || 
						    minor_count == SEQMAXUNSIGNED))
		      ? Counts.g : minor_count;
                    minor_count = (Counts.c > 0 && (Counts.c < minor_count || 
						    minor_count == SEQMAXUNSIGNED))
		      ? Counts.c : minor_count;
                    minor_count = (Counts.t > 0 && (Counts.t < minor_count || 
						    minor_count == SEQMAXUNSIGNED))
		      ? Counts.t : minor_count;
                    minor_count = (Counts.zero > 0 && (Counts.zero < minor_count || 
						       minor_count == SEQMAXUNSIGNED))
		      ? Counts.zero : minor_count;
                    minor_count = (Counts.one > 0 && (Counts.one < minor_count || 
						      minor_count == SEQMAXUNSIGNED))
		      ? Counts.one : minor_count;
                  }
                else if (freqfilter == 0 )
                  minor_count = 1;

                if(strict == false && numstates > 1 && minor_count > freqfilter)
                  //site is variable
                  {
                    _positions.push_back(double(j+1));//add 1 to emulate a real positions
                  }
                else if(strict == true && numstates == 2 && minor_count > freqfilter )
                  //if there are 2 and only 2 character states
                  {
                    _positions.push_back(double(j+1));//add 1 to emulate a real positions
                  }
              }
          }
      }
    std::vector<std::string> _data(numseqs);
    for(unsigned i = 0 ; i < numseqs ; ++i)
      {
        for(unsigned j = 0 ; j < _positions.size() ; ++j)
          {
            if (alignment[i][unsigned(_positions[j])-1] == IDENTICAL)
              {
                _data[i] += alignment[0][unsigned(_positions[j])-1];
              }
            else
              {
                _data[i] += alignment[i][unsigned(_positions[j])-1];
              }
          }
      }
    PolyTable::assign(&_positions[0],_positions.size(),&_data[0],_data.size());
  }
}
