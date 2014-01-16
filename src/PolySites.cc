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

#include <Sequence/PolySites.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/Portability/StringStreams.hpp>
#include <iterator>
#include <algorithm>
#include <iostream>

/*! \class Sequence::PolySites Sequence/PolySites.hpp
  \ingroup polytables
  This is one of the more useful classes in namespace Sequence.
  Its purpose is to take a bunch of data (a vector<Fasta> in fact),
  and turn it into a list of variable positions.  It doesn't matter
  whether or not you have an outgroup in you vector (except that if you
  want to use it for later analysis, it had better be present).\n
  \n
  The default behavior of this class is just to play with the std::strings
  themselves.  So what you end up with is a vector of variable sites,
  stored in PolyData::positions, and a vector ofstd::strings containing 
  the variable sites, stored in PolyData::data.  Note that if you include
  an outgroup in your vector<Fasta>, and it contains a different character
  than the ingroup at some site, then the site is considered variable.\n
  \n
  You can also try and turn the data into a "binary" (i.e. 0 and 1) format, by
  a call to Sequence::PolySites::Binary.
  \n
  EXAMPLE:\n
  Here is a common use of a PolySites class.  You have a file, "gene.fasta", containing
  some number of sequences that represent polymorphism data.  The file is assumed to be
  aligned, but we'll check for that, just in case you forgot to run ClustalW or something.\n
  \n
  \code
  #include <string>
  #include <iostream>
  #include <Sequence/Fasta.hpp>
  #include <Sequence/Alignment.hpp>
  #include <Sequence/SeqExceptions.hpp>
  int main(int argc, char *argv[]) {
  const char *infile = "gene.fasta";
  vector<Fasta> data;
  try {
  Sequence::Alignment::GetData (data,infile);
  assert( Sequence::Alignment::IsAlignment (data) );
  if ( Sequence::Alignment::Gapped (data) ) 
  Sequence::Alignment::RemoveTerminalGaps (data);
  }
  catch (SeqException &e)
  {
  cerr << "uh-oh! processing file gene.fasta resulted in throwing an exception"<<endl;
  e.print(cerr);
  cerr << endl;
  exit(1);
  }
  PolySites *polytable = new PolySites(data);
  }
  \endcode
  \n
  Removing the terminal gaps guarantees that polymorphic site positions are labelled
  starting from the first ungapped position. Of course, a lot of the extra syntax
  in the example can be eliminated by giving the following 2 using declarations:\n
  using namespace Sequence;\n
  using namespace Sequence::Alignment;\n
  \n
  For a second example, assume the data are in the file "gene.aln", the results of a
  ClustalW alignment.\n
  \n
  \code
  #include <iostream>
  #include <Sequence/Clustalw.hpp>
  #include <Sequence/Fasta.hpp>
  #include <Sequence/SeqExceptions.hpp>
      
  int main(int argc, char *argv[]) {
  istream in;
  in.open("gene.aln");
  ClustalW<Fasta> aligned_data;
  try {
  in >> aligned_data;
  assert(aligned_data.IsAlignment());
  if(aligned_data.Gapped())
  aligned_data.RemoveTerminalGaps();
  }
  catch (SeqException &e)
  {
  cerr << "uh-oh! processing file gene.aln resulted in throwing an exception"<<endl;
  e.print(cerr);
  cerr << endl;
  exit(1);
  }
  PolySites *polytable = new PolySites(aligned_data.GetData());
  }
  \endcode
  @short Polymorphism tables for sequence data
*/

namespace Sequence
{
  PolySites::PolySites (void)
  {}

  PolySites::PolySites (const std::vector < double > &List, const std::vector <std::string > &stringList):
    PolyTable()
    /*!
      Use this constructor if you already have a list of positions and characters
      \param List a list of doubles representing positions of polymorphic positions
      \param stringList a vector of strings representing the polymorphic characters
    */
  {
    seqlen = stringList[0].length();
    numseqs = stringList.size();
    PolyTable::assign(&List[0],List.size(),&stringList[0],stringList.size());
  }

  PolySites::PolySites (PolyTable::const_site_iterator beg,
			PolyTable::const_site_iterator end):
    PolyTable(beg,end)
  {
  }

  std::istream & PolySites::read(std::istream &s) 
    

  {
    std::vector<double> _pos;
    std::vector<std::string> _ind;
    std::string temp;
    //get positions--which is the first line of input
    std::getline(s,temp);
    istr i(temp);
    std::copy(std::istream_iterator<double>(i),std::istream_iterator<double>(),
	      std::back_inserter(_pos));

    while (std::getline(s,temp))
      {
	std::string temp2;
	istr i(temp);
	std::copy(std::istream_iterator<char>(i),std::istream_iterator<char>(),
		  std::back_inserter(temp2));
	_ind.push_back(temp2);
      }
    if (! assign(&_pos[0],_pos.size(),&_ind[0],_ind.size()) )
      {
	throw badFormat("PolySites::read() -- format error, unable to assign data");
      }
    return s;
  }

  std::ostream & PolySites::print(std::ostream &stream) const
  /*!
    Allows objects of type Sequence::PolySites to be 
    written to output streams.  The output is a simple,
    tab-delimited table of variable site positions and
    characters
    \note segsite positions are output with the count starting from 1, not zero
    @brief output a tab-delimited array of positions and character states
  */
  {
    for(unsigned i = 0 ; i < this->numsites() ; ++i)
      {
        if(i==0)
          stream << this->position(i);
        else
          stream << '\t' << this->position(i);
      }
    stream << '\n';

    for(unsigned i = 0 ; i < this->size() ; ++i)
      {
        for(unsigned j = 0 ; j < this->numsites() ; ++j)
          {
            if(j==0)
              stream << (*this)[i][j];
            else
              stream << '\t'<<(*this)[i][j];
          }
        if(i < this->size()-1)
          stream << '\n';
      }
    return (stream);
  }
}
