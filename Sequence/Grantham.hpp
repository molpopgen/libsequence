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

#ifndef GRANTHAM_H
#define GRANTHAM_H
/*! \file Grantham.hpp
  @brief Grantham's distances (Sequence::Grantham)
*/

/*! \class Sequence::Grantham Sequence/Grantham.hpp
  A functor to return the Grantham's distance between
  two amino acids.
 
  @short Grantham's distances
*/

namespace Sequence
  {
  class Grantham
    {
    private:
      double D[60][60];
    public:
      Grantham(void);
      double operator()(char aa1, char aa2) const;
    };
}
#endif
