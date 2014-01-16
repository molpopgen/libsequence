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

#ifndef __SEQ_PROPERTIES_HPP__
#define __SEQ_PROPERTIES_HPP__
/*! \file SeqProperties.hpp
  \brief Function objects to ask simple questions about sequence composition.
  Declares Sequence::ambiguousNucleotide,Sequence::invalidPolyChar.
 */
#include <functional>
#include <cctype>
#include <Sequence/Comparisons.hpp>

namespace Sequence
{
  struct ambiguousNucleotide : public std::unary_function<char,bool>
			       /*! \struct ambiguousNucleotide Sequence/SeqProperties.hpp
				*/
  {
    inline bool operator()(const char & c) const
    /*!
      \return true if c is not A,G,C, or T, false otherwise
      \note Case-insensitive
    */
    {
      const char ch = char(std::toupper(c));
      return (ch != 'A' &&
	      ch != 'G' &&
	      ch != 'T' &&
	      ch != 'C' );
    }
  };

  struct invalidPolyChar : public std::unary_function<char,bool>
			   /*! \struct invalidPolyChar Sequence/SeqProperties.hpp
			     This functor can be used to determine
			     if a range contains characters that
			     the SNP analysis routines in this
			     library cannot handle gracefully
			    */
  {
    inline bool operator()(const char & nucleotide) const
    /*!
      \return true if c is not in the set {A,G,C,T,N,-,.}, false otherwise. The period
      (.) can be used as an "identical to the 1st seq in a file" character, so 
      should be considered valid
      \note Case-insensitive
    */
    {
      const char ch = char(std::toupper(nucleotide));
      return ( ch != 'A' &&
	       ch != 'G' &&
	       ch != 'C' &&
	       ch != 'T' &&
	       ch != 'N' &&
	       ch != '-' && 
	       ch != '.' );
    }
  };

}

#endif
