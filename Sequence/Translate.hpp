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

#ifndef __TRANSLATE_HPP__
#define __TRANSLATE_HPP__
#include <string>
#include <Sequence/SeqEnums.hpp>
#include <Sequence/SeqExceptions.hpp>
/*! \file Translate.hpp
  @brief declares Sequence::Translate,a function to translate CDS sequences into peptide sequences
*/

/*!
  \defgroup misc Miscellany
 */
namespace Sequence
  {
  /*!
    \ingroup misc
    \param beg a pointer to the beginning of the region to translate
    \param end a pointer to 1 past the end of the region to translate
    \param genetic_code must be a value from the enumeration list Sequence::GeneticCodes
    \param gapchar a character representing an alignment gap
    \return a string representing the translation of the range
    \throw Sequence::SeqException if \a genetic_code is invalid
    \code
    #include <Sequence/Translate.hpp>
    \endcode
  */
    std::string Translate(std::string::const_iterator beg,
			  std::string::const_iterator end,
			  Sequence::GeneticCodes 
			  genetic_code = Sequence::UNIVERSAL,
			  const char & gapchar = '-')
      ;
}
#endif
