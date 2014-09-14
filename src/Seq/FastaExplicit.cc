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

#include <Sequence/Fasta.hpp>
#include <Sequence/Clustalw.hpp>
#include <Sequence/phylipData.hpp>
#include <Sequence/Alignment.hpp>
/*! \file FastaExplicit.cc
  For details, see documentation of FastaExplicit.hpp
  @brief Explicit instantiations of templates for type Sequence::Fasta
*/
namespace Sequence
{
  /*! An explicit instantation for type Sequence::Fasta for GNU systems.
    See file documentation for FastaExplicit.hpp
  */
  template class Sequence::ClustalW< Sequence::Fasta >;
  /*! An explicit instantation for type Sequence::Fasta for GNU systems.
    See file documentation for FastaExplicit.hpp
  */
  template class Sequence::phylipData< Sequence::Fasta >;
  namespace Alignment
  {
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template void GetData( std::vector<Sequence::Fasta> &, const char *);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template std::istream & GetData( std::vector<Sequence::Fasta> &,
				     std::istream &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template std::istream & ReadNObjects( std::vector<Sequence::Fasta> &, 
					  unsigned, std::istream &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template bool Gapped(const std::vector<Sequence::Fasta> &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template bool IsAlignment(const std::vector<Sequence::Fasta> &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template bool validForPolyAnalysis(std::vector<Sequence::Fasta>::const_iterator beg,
				       std::vector<Sequence::Fasta>::const_iterator end);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template unsigned UnGappedLength(const std::vector<Sequence::Fasta> &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template void RemoveGaps(std::vector<Sequence::Fasta> &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template void RemoveTerminalGaps(std::vector<Sequence::Fasta> &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template void RemoveFixedOutgroupInsertions( std::vector<Sequence::Fasta> & data,
						 unsigned site,
						 const unsigned & ref );
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template std::vector<Sequence::Fasta>
    Trim(const std::vector<Sequence::Fasta> &,
	 const std::vector<int> &);
    /*! An explicit instantation for type Sequence::Fasta for GNU systems.
      See file documentation for FastaExplicit.hpp
    */
    template std::vector<Sequence::Fasta>
    TrimComplement(const std::vector<Sequence::Fasta> &,
		   const std::vector<int> &);
  }

}
