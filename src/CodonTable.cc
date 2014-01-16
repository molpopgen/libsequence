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

/*! \example codons.cc
 */
#include <Sequence/Seq.hpp>
#include <Sequence/CodonTable.hpp>
#include <Sequence/Comparisons.hpp>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iostream>
#include <cctype>
#include <vector>
#include <string>
#include <utility>

using ::toupper;

namespace
{
  typedef std::pair< std::string, int > CodonFreq;
  //typedef std::string::size_type sst;

  //our alphabet consists of nucleotide characters
  static const unsigned int alphsize = 4;
  static const char alphabet[alphsize] ={'A','C','G','T'};

  //Sequence::CodonUsageTable doWork(const std::string &sequence)
  Sequence::CodonUsageTable doWork( std::string::const_iterator beg,
				    std::string::const_iterator end )
  {
    Sequence::CodonUsageTable UsageTable;
    for (unsigned i = 0 ; i < alphsize ; ++i)
      for (unsigned j = 0 ; j < alphsize ; ++j)
        for (unsigned k = 0 ; k < alphsize ; ++k)
          {
            std::string codon;
            codon += alphabet[i];
            codon += alphabet[j];
            codon += alphabet[k];
            size_t pos = 0;
            unsigned count = 0;
	    for( ; (beg+pos) < end ; pos += 3 )
	      {
		//Ask if positions beg to beg+3 differ from codon, but do 
		//not allow missing data to result in the answer to be false
		if( ! Sequence::Different( std::string(beg+pos,beg+pos+3),codon,false) )
		  {
		    ++count;
		  }
	      }
            UsageTable.push_back( CodonFreq(codon,count) );
          }
    return UsageTable;
  }
}  

namespace Sequence
{
  CodonUsageTable makeCodonUsageTable(const Seq* sequence)
  /*!
    \param sequence and object of type Sequence::Seq1
    \return and object of type Sequence::CodonUsageTable
    \note Assumes first character of sequence is a first codon position.
  */
  {
    return doWork(sequence->begin(),sequence->end());
  }

  CodonUsageTable makeCodonUsageTable(const std::string &sequence)
  /*!
    \param sequence and object of type std::string
    \return and object of type Sequence::CodonUsageTable
    \note Assumes first character of sequence is a first codon position
  */
  {
    return doWork(sequence.begin(),sequence.end());
  }

  CodonUsageTable makeCodonUsageTable(std::string::const_iterator beg,
				      std::string::const_iterator end)
  /*!
    \param beg a const_iterator to the beginning of a std::string or Sequence::Seq
    \param end a const_iterator to the end of a std::string or Sequence::Seq
    \return and object of type Sequence::CodonUsageTable
    \note beg and end can be adjusted to point to the first at last positions in a CDS
  */
  {
    //we do the typecast to string because it allows the use of std::string::find,
    //which is efficient for searching
    return doWork(beg,end);
  }
}
