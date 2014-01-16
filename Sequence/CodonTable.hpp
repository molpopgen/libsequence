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

#ifndef _CODON_TABLE_H_
#define _CODON_TABLE_H_
#include <Sequence/typedefs.hpp>

/*! \file CodonTable.hpp
  @brief facility to count codons in CDS sequence, function Sequence::makeCodonUsageTable
*/

namespace Sequence
{
  /*!
    \c #include \c <Sequence/CodonTable.hpp>
    A codon usage table is a list of the codons and
    the number of times each codon occurs in a sequence.
    This is represented by the type Sequence::CodonUsageTable,
    which is a vector < pair<string,int >.
    To output a codon table:
    \code
    Sequence::Fasta sequence;
    cin >> sequence;
    CodonUsageTable x = makeCodonUsageTable(&sequence);
    for(unsigned i = 0 ; i < x.size ; ++i)
    {
    cout << x[i].first << '\t' << x[i].second << '\n';
    }
    \endcode
  */
  class Seq;
  CodonUsageTable makeCodonUsageTable(const Seq* sequence);
  CodonUsageTable makeCodonUsageTable(const std::string &sequence);
  CodonUsageTable makeCodonUsageTable(std::string::const_iterator beg,
				      std::string::const_iterator end);

}
#endif
