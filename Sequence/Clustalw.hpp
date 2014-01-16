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

#include <Sequence/AlignStream.hpp>

#ifndef __CLUSTALW_HPP__
#define __CLUSTALW_HPP__

/*! \file Clustalw.hpp
  @brief Declaration of Sequence::ClustalW streams
*/

/*! \class Sequence::ClustalW Sequence/Clustalw.hpp
  \ingroup alignment
  This class defines an input routine for alignments in ClustalW
  format.
  This class template can be instantiated with the following types:
  \code
  Sequence::ClustalW< Sequence::Fasta >; //or anything derived from Sequence::Seq
  Sequence::ClustalW< std::pair<std::string, std::string> >;
  \endcode
  \code
  //example: a std::pair<std::string,std::string>
  #include<utility>
  #include <string>
  #include <iostream>
  #include <Sequence/Clustalw.hpp>
  
  //using namespace std;
  using namespace Sequence;
  
  int main()
  {
  //instantiate the template
  ClustalW< std::pair<string,string> > x;
  //read in an object from stdin
  cin >> x;
  for (int i = 0 ; i < x.size() ; ++i)
  {
  //output each element in Fasta format
  //but without using Sequence::Fasta
  cout << '>'
  << x[i].first
  << '\n'
  << x[i].second
  << '\n';
  }
  }
  \endcode
  
  @short ClustalW streams
*/

namespace Sequence
{
  template < typename T >
  class ClustalW:public AlignStream < T >
  {
  private:
    std::istream & ReadThroughLine (std::istream &);
  public:
    ClustalW (void):AlignStream<T>()
    {}
    ClustalW(const AlignStream<T> & a) : AlignStream<T>(a)
    {
    }
    ~ClustalW(void)
    {}
    std::istream & read (std::istream & s);
    std::ostream & print(std::ostream &s) const;
    ClustalW<T> & operator=( const AlignStream<T> & rhs)
    {
      this->assign(rhs.begin(),rhs.end());
      return *this;
    }
  };
}
#include <Sequence/bits/Clustalw.tcc>
#endif
