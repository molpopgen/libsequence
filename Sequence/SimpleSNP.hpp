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

#ifndef __SIMPLE_SNP_HPP__
#define __SIMPLE_SNP_HPP__
/*! \file SimpleSNP.hpp
  @brief  Declaration of Sequence::SimpleSNP, a polymorphism table stream in a "spreadsheet" format
*/
/*! \class Sequence::SimpleSNP Sequence/SimpleSNP.hpp
  \ingroup polytables
  This class is used to deal with simple "polymorphism tables" that are
  formatted in a particular way.  

  The major purpose of this class is to be able to read in data files
  in this format.  An example of the format is as follows:\n
  \n
  6	6\n
  679	1004	1153	1155	1277	1295	\n
  N	N	N	N	N	N	\n
  Sim2	A	C	C	G	A	A\n
  Sim3	G	C	T	G	A	C\n
  Sim4	A	C	C	T	G	C\n
  Sim5	A	C	C	T	G	C\n
  Sim7	A	T	C	T	A	C\n
  Sim8	A	C	C	G	A	A\n
  \n
  The two numbers on the first line are the sample size and number of 
  segregating sites, respectively.  The next line contains the positions
  of each site.  The third line contains the outgroup state for each variable
  site--if this is unknown, use an 'N' to indicate the ambiguity. 
  The rest of the lines contain a unique
  sequence name, and then the states of each segregating site. 
 
  Finally, if all of the outgroup states are missing (as is the case in the 
  above example, then it is assumed that no outgroup was typed, and no outgroup 
  sequence is stored in PolyTable::data.  However, if at least one of the outgroup
  states is unambiguous, a std::string representing outgroup states is stored as the first
  element in PolyTable::data.
 
  @short SNP table data format
*/

#include <Sequence/PolyTable.hpp>

namespace Sequence
  {
  class SimpleSNP:public PolyTable
    {
    private:
      mutable bool Diploid,isoFemale;
      bool haveOutgroup;
      std::vector<std::string> _names;
    public:
      SimpleSNP (const bool diploid =0,const bool isofemale=0)  : PolyTable(),
          Diploid(diploid),isoFemale(isofemale),haveOutgroup(false)
          /*!
          The two bools that this constructor takes allow you to deal
          with two very different types of polymorphism data.  If both
          bools are set to 0 (the default), the data are simply read in
          as they are.  However, if diploid == 1, then if n sequences are
          read in, the data are assumed to be diploid (make sense...), and
          are converted into 2n strings.  Further, if heterozygous bases
          are encoutered (R,W, etc.), the two possible states are arbitrarily
          assigned to each sequence.\n
          \n
          If isofemale==1, it is assumed that the
          data represent real haplotype data (i.e. phase is known).  The name
          for the bool comes from the fact that data gathered from Drosophila 
          lines is often obtained from isofemale stocks, making them homozygous
          such that the phase of each SNP is known.  If a heterozygous base is
          found, one of the two possible states will be assigned randomly
          (NOT IMPLEMENTED YET!)
          */
      {}
      ~ SimpleSNP (void) {}
     
      bool outgroup(void) const;
      void set_outgroup( const bool & b );
      std::string label(unsigned i) const;
      std::istream & read (std::istream & s) ;
      std::ostream & print(std::ostream &o) const;
    };
    typedef SimpleSNP Hudson2001;
}
#endif
