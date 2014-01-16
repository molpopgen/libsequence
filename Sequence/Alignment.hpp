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

#ifndef __ALIGNMENT_HPP__
#define __ALIGNMENT_HPP__
#include <iosfwd>
#include <fstream>
#include <vector>
#include <string>
#include <Sequence/SeqExceptions.hpp>

/*! \file Alignment.hpp
  @brief Declaration of namespace Sequence::Alignment
*/

/*!
  \defgroup alignment Classes and functions dealing with alignments
  \ingroup data
*/
/*! \namespace Sequence::Alignment
  \ingroup alignment
  \#include<Sequence/Alignment.hpp>
  \n
  This namespace defines a series of template routines used for manipulating aligned
  sequence data. The routines may be used on their own, and are also the building 
  blocks for alignment I/O as defined in Sequence::AlignStream.

  \note Valid types for instantiation of these templates are the following:
  std::pair<std::string,std::string>, or anything derived from there, which 
  include Sequence::Seq, etc.  This requirment is checked at compile time 
  via compile-time assertions using the <A HREF="http://www.boost.org">boost library</A>.
  There are also specializations such that std::string is also a valid type.
  @short Routines fundamental to aligned data
*/
namespace Sequence
{
  namespace Alignment
  {
    //fill a vector<T> from a file
    template < typename T >
    void GetData (std::vector < T >&seqarray,
                  const char *infilename);
    //fill a vector<T> from a stream
    template < typename T >
    std::istream & GetData (std::vector < T >&seqarray, 
			    std::istream & input_stream);
    //fill a vector<T> with N objects from a stream
    template < typename T >
    std::istream & ReadNObjects ( std::vector < T > &seqarray,   unsigned n,
				  std::istream & input_stream);
    //delete a vector of pointers
    template < typename T >
    void EmptyVector (std::vector< T * > &seqarray);
    //check if elements of data contain a gap character '-'
    template < typename T >
    bool Gapped (const std::vector < T >&data);
    //verify that all elements of data are the same length
    template < typename T >
    bool IsAlignment (const std::vector < T  >&data);
    //check to make sure that an aligment contains 
    //characters that the SNP routines in this library can handle
    template<typename Iterator>
    bool validForPolyAnalysis( Iterator beg,
			       Iterator end );
    //return the length of the alignment without gaps
    template < typename T >
    unsigned UnGappedLength (const std::vector <T>&data) ;
    //remove all gaps from an alignment
    template <typename T>
    void RemoveGaps (std::vector <T> &data);
    //remove all gaps from the ends of an alignment
    template < typename T >
    void RemoveTerminalGaps (std::vector < T >&data);
    template < typename T >
    void RemoveFixedOutgroupInsertions( std::vector<T> & data,
				unsigned site,
				const unsigned & ref );
    //see manual
    template < typename T >
    std::vector < T >Trim (const std::vector < T >&data,
                           const std::vector <int> &sites) 
      ;
    //see manual
    template < typename T >
    std::vector < T >TrimComplement (const std::vector < T >&data,
                                     const std::vector <int> &sites) 
      ;

    //declaration of specializations
    template<> bool Gapped(const std::vector<std::string> &data);
    template <> bool IsAlignment(const std::vector<std::string> &data);
    template<>
    bool validForPolyAnalysis( std::vector<std::string>::const_iterator beg,
 			       std::vector<std::string>::const_iterator end );
    template<>
    bool validForPolyAnalysis( std::vector<std::string>::iterator beg,
			       std::vector<std::string>::iterator end );
    template <> 
    unsigned UnGappedLength (const std::vector <std::string>&data) ;
    template <>
    void RemoveGaps (std::vector <std::string> &data);
    template <>
    void RemoveTerminalGaps (std::vector <std::string> &data);
    template <>
    void RemoveFixedOutgroupInsertions(std::vector<std::string> &data,
				       unsigned site,
				       const unsigned & ref);
    template <>
    std::vector <std::string >Trim (const std::vector <std::string>&data,
                                    const std::vector <int> &sites) ;
    
    template <>
    std::vector <std::string> TrimComplement (const std::vector <std::string>&data,
					      const std::vector <int> &sites) ;

  }
}
#include <Sequence/bits/Alignment.tcc>
#endif
