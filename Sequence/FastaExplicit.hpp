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

#ifndef __FASTA_EXPLICIT_HPP_
#define __FASTA_EXPLICIT_HPP_

/*! \file FastaExplicit.hpp
  This file contained explicit instantiations for functions
  in namespace Sequence::Alignment and for class Sequence::ClustalW.
  The type of the instantiation is Sequence::Fasta.\n
  Example of how to use this portably:
  \code
  //suppose you want to use features from
  //both namespace Sequence::Alignment and
  //use the Sequence::ClustalW template
 
  //check for GNU C++ of at least version 3
  #if defined( __GNUG__ ) && __GNUC__ >= 3
  //take advantage of the explicit instantations
  //at link-time
  #include <Sequence/FastaExplicit.hpp>
  #else
  //need to do implicit instanstiation
  #include <Sequence/Fasta.hpp>
  #include <Sequence/Alignment.hpp>
  #include <Sequence/Clustalw.hpp>
  #endif
  \endcode
  \warning The instantiations make use of a GNU C++ extension
  to the rules for explicit instantiation of templates. The declarations
  in this file are "extern", and their definitions are in FastaExplicit.cc.
  These declarations may not be portable. In fact, compiling with -ansi
  and -pendantic will result in a compiler error.  To prevent the error:
  \code
  #if defined( __GNUG__ ) && __GNUC__ >= 3 && (!__STRICT_ANSI__)
  #include <Sequence/FastaExplicit.hpp>
  #else
  #include <Sequence/Fasta.hpp>
  #include <Sequence/Alignment.hpp>
  #include <Sequence/Clustalw.hpp>
  #endif
  \endcode
  @brief delcaration of explicit instantiations of namespace Sequence::Alignment
  and Sequence::ClustalW for type Sequence::Fasta
*/
//include the relevant headers
#if __GNUG__ && __GNUC__ >= 3
#include <Sequence/Clustalw.hpp>
#include <Sequence/phylipData.hpp>
#include <Sequence/Fasta.hpp>
namespace Sequence
{
  namespace Alignment
  {
    extern template void GetData( std::vector<Sequence::Fasta> &,
				  const char *);
    extern template std::istream & GetData( std::vector<Sequence::Fasta> &,
					    std::istream & );
    extern template std::istream & ReadNObjects( std::vector<Sequence::Fasta> &,
						 unsigned, std::istream &);
    extern template bool Gapped(const std::vector<Sequence::Fasta> &);
    extern template bool IsAlignment(const std::vector<Sequence::Fasta> &);
    extern template bool validForPolyAnalysis(std::vector<Sequence::Fasta>::const_iterator beg,
					      std::vector<Sequence::Fasta>::const_iterator end);
    extern template unsigned UnGappedLength(const std::vector<Sequence::Fasta> &);
    extern template void RemoveGaps(std::vector<Sequence::Fasta> &);
    extern template void RemoveTerminalGaps(std::vector<Sequence::Fasta> &);
    extern template void RemoveFixedOutgroupInsertions( std::vector<Sequence::Fasta> & data,
							unsigned site,
							const unsigned & ref );
    extern template std::vector<Sequence::Fasta>
    Trim(const std::vector<Sequence::Fasta> &, const std::vector<int> &);
    extern template std::vector<Sequence::Fasta>
    TrimComplement(const std::vector<Sequence::Fasta> &, const std::vector<int> & );
  }
  extern template class Sequence::ClustalW< Sequence::Fasta >;
  extern template class Sequence::phylipData< Sequence::Fasta >;
}
#endif

#endif
