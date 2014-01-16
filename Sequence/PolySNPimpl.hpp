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
#include <Sequence/stateCounter.hpp>
#include <string>

namespace Sequence
{
  struct _PolySNPImpl
  /*!
    Implementation details for PolySNP.  This class is visible
    so that it can be accessed from classes derived from PolySNP.
    A PolySNP object contains a pointer to an instance of this class
    that is storage class protected.
  */
  {
    const PolyTable* _data;
    unsigned _nsites,_nsam,_outgroup;
    bool _haveOutgroup, _totMuts;
    unsigned _totsam;
    unsigned _DVK;
    double _DVH;
    bool _counted_singletons;
    bool _know_pi;
    bool _CalculatedDandV;
    double _pi;
    unsigned _singletons;
    unsigned _walls_Bprime,_NumPoly;
    double _walls_B,_walls_Q;
    bool _calculated_wall_stats;
    std::vector< Sequence::stateCounter > _counts;
    std::vector< std::pair< bool, Sequence::stateCounter > > _derivedCounts;
    bool _preprocessed;
    void preprocess(void);

    _PolySNPImpl (const Sequence::PolyTable * data, bool haveOutgroup ,
		  unsigned outgroup, bool totMuts);
  };
}
