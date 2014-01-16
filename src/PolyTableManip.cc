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

#include <Sequence/PolyTable.hpp>
#include <Sequence/PolyTableManip.hpp>

namespace Sequence
{
  polySiteVector rotatePolyTable(const Sequence::PolyTable *data)
  /*!
    Rotate a polymorphism table
    into a vector of pairs, where the
    pairs are of type std::pair<double, string>,
    representing the site position and the characters
    at that site
    \param data a pointer to a Sequence::PolyTable
    \ingroup polytables 
  */
  {
    polySiteVector L;
    for (unsigned i = 0 ; i < data->numsites() ; ++i)
      {
        std::string s;
        for(unsigned j = 0 ; j < data->size() ; ++j)
          {
            s += (*data)[j][i];
          }
        L.push_back( polymorphicSite(data->position(i), s) );
      }
    return L;
  }
}
