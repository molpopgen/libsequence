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

#ifndef SINGLESUB_H
#define SINGLESUB_H
/*! \file SingleSub.hpp
  @brief used by Sequence::Comeron95, class Sequence::SingleSub calculates divergence between codons that differ at one site
*/

/*!
  \class Sequence::SingleSub Sequence/SingleSub.hpp
  \ingroup kaks
  A functor to obtain divergence statistics for Comeron's method for codons that differ at one position.  Used by
  Sequence::Comeron95

  @author Kevin Thornton
  @short Deal with codons differing at 1 position
*/
#include <string>

namespace Sequence
  {
  class RedundancyCom95;
  class SingleSub
    {
    private:
      double q0i, q2Si, q2Vi, q4i, q0j, q2Sj, q2Vj, q4j, p0i, p2Si, p2Vi,
      p4i, p0j, p2Sj, p2Vj, p4j;
      void Calculate (const RedundancyCom95 * sitesObj, const std::string & cod1,
                      const std::string & cod2);
    public:
      explicit SingleSub(void)
      {}
      void operator()(const RedundancyCom95 * sitesObj,
                      const std::string &cod1,
                      const std::string &cod2);
      ~SingleSub (void)
      {}
      double P0(void) const;
      double P2S(void) const;
      double P2V(void) const;
      double P4(void) const;
      double Q0(void) const;
      double Q2S(void) const;
      double Q2V(void) const;
      double Q4(void) const;
  };
}
#endif
