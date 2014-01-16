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

#include <Sequence/SeqExceptions.hpp>
#include <iostream>

namespace Sequence
{
  SeqException::SeqException (void):message("Generic SeqException: error encountered")
    /*!
      Default constructor--will generate a default error message,
      which is not likely to be meaningful
      \note this constructor only exists for the purpose of inheriting from this class
    */
  {
  }
  SeqException::SeqException (const char *x):message(x)
    /*!
      Throw the exception with error message x
    */
  {}
  SeqException::~SeqException(void)
  {}
  std::ostream & SeqException::print (std::ostream & out)
    /*!
      Write the error to out.
    */
  {
    out << message;
    return out;
  }
  std::ostream & SeqException::print (std::ostream & out) const
    /*!
      Write the error to out.
    */
  {
    out << message;
    return out;
  }
  const char* SeqException::error (void) const
    /*!
      Return the error message
      if you want to use it in
      some other fashion then
      printing it to an std::ostream
    */
  {
    return message;
  }
    
  badFormat::~badFormat(void)
  {}
  badFormat::badFormat (const char * x):SeqException(x)
    /*!
      Throw the exception with error message x
    */
  {
    return;
  }
}
