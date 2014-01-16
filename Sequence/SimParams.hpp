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

#ifndef SIMPARAMS_H
#define SIMPARAMS_H
/*! \file SimParams.hpp
  @brief Sequence::SimParams reads in the parameters of Dick Hudon's coalescent simulation program.  Used in conjunction with Sequence::SimData
*/

/*! \class Sequence::SimParams Sequence/SimParams.hpp
  \ingroup coalescent
  include SimParams.h
  Allows reading in and printing out of the parameter
  list that Hudson's coalescent simulation program spits
  out at the beginning of its execution.  An example of use
  is found in tajd.cc in the Examples section.
 
  @author Kevin Thornton
  @short Parameters for Hudson's simulation program
*/
#include <iosfwd>
#include <string>
#include <vector>
#include <cstdio>

namespace Sequence
{
  class SimParams
  {
    friend std::ostream& operator<<(std::ostream&,class SimParams &object);
  private:
    std::string _command_line;
    unsigned _howmany, tsam;
  public:
    SimParams(void);
    std::istream& read(std::istream& s);
    int fromfile ( FILE * openfile );
    std::string params (void) const
      /*!
	\return the command-line input to ms
	\note for complicated models, this can be parsed
	with a stringstream to figure out what the parameters are
      */
    {
      return _command_line;
    }
    unsigned totsam (void) const
      /*!
	\return the total sample size (# gametes)
      */
    {
      return (tsam);
    }
    unsigned runs (void) const
      /*!
	\return number of genealogies to generate
      */
    {
      return (_howmany);
    }
  };

  std::istream& operator>>(std::istream& s,  SimParams& c);
}
#endif
