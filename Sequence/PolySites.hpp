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

#ifndef POLYSITES_H
#define POLYSITES_H
/*! \file PolySites.hpp
  @brief Sequence::PolySites, generates polymorphism tables from data
*/
#include <Sequence/PolyTable.hpp>
namespace Sequence
  {
  class Fasta;
  class PolySites : public PolyTable
    {
    private:
      /*!
      PolySites::fillIt() is the function that actually fills the polymorphism table.
      */
      template<class __DataType>
      void fillIt(const std::vector < __DataType >&alignment,
                         bool strictInfSites = 0,
                         bool ignoregaps = 1,bool skipMissing=false,
                         unsigned freqfilter=0);
    protected:
      size_t numseqs;
      size_t seqlen;
    public:
      PolySites (void);
      template<typename __DataType>
      PolySites (const std::vector < __DataType >&alignment, bool strictInfSites =
                   0, bool ignoregaps = 1,bool skipMissing=false,
                 bool skipAdjSNP=false, unsigned freqfilter=0);
      PolySites (const std::vector < double > &List, const std::vector < std::string > &stringList);
      PolySites (PolyTable::const_site_iterator beg,
		 PolyTable::const_site_iterator end);
      ~PolySites(void)
      {}
      std::istream & read(std::istream &s) ;
      std::ostream & print(std::ostream &stream) const;
    };
}
#include <Sequence/bits/PolySites.tcc>
#endif
