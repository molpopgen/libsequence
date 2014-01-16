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

#include <Sequence/SeqConstants.hpp>
#if ( defined(HAVE_CLIMITS) && (!(defined(HAVE_LIMITS))) )
#include <climits>
#else
#include <limits>
#endif

namespace Sequence
{
#if ( defined(HAVE_CLIMITS) && (!(defined(HAVE_LIMITS))) )
  /*! \var const unsigned SEQMAXUNSIGNED
    The maximum value of an unsinged integer.
  */
  const unsigned SEQMAXUNSIGNED = UINT_MAX;
  /*! \var const unsigned SEQMAXDOUBLE
    The maximum value of an double
  */
  const double SEQMAXDOUBLE = DOUBLE_MAX;
#else
  /*! \var const unsigned SEQMAXUNSIGNED
    The maximum value of an unsinged integer.
  */
  const unsigned SEQMAXUNSIGNED = std::numeric_limits<unsigned>::max();
  /*! \var const unsigned SEQMAXDOUBLE
    The maximum value of an double
  */
  const double SEQMAXDOUBLE = std::numeric_limits<double>::max();
#endif
}
