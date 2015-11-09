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

#ifndef __FST_H__
#define __FST_H__

#include <utility>
#include <set>
#include <memory>
#include <Sequence/stateCounter.hpp>
#include <Sequence/SeqExceptions.hpp>

/*! \file FST.hpp
  @brief delcaration of a class (Sequence::FST) to analyze population structure

  \deprecated Will be removed in libsequence 2.0
*/
namespace Sequence
{
  class PolyTable;
  struct FSTimpl;
  class FST 
  {
  private:
    void doCalcs(void) const;
    std::unique_ptr<FSTimpl> impl;
  public:
    explicit FST(const PolyTable *data, unsigned npop, const unsigned *config=NULL,
		 const double *weights=NULL, bool haveOutgroup = false,
		 unsigned outgroup = 0);
    FST(const FST &) = delete;
    FST & operator=(const FST &) = delete;
    ~FST(void);
    double HSM(void) const;
    double Slatkin(void) const;
    double HBK(void) const;
    double piB(void) const;
    double piT(void) const;
    double piS(void) const;
    double piD(void) const;
    std::set<double> shared(unsigned pop1, unsigned pop2) const;
    std::set<double> fixed(unsigned pop1, unsigned pop2) const;
    std::pair< std::set<double>,std::set<double> > Private(unsigned pop1, unsigned pop2) const;
  };
}
#endif
