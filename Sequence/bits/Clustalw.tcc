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

// Code for the -*- C++ -*- namespace Sequence::ClustalW<T>

/*! \file Clustalw.tcc
  @brief code for Clustalw.hpp
*/
#include <map>
#include <Sequence/AlignStream.hpp>
#include <iterator>
#include <algorithm>

namespace Sequence
{
  template < typename T >
  std::istream & ClustalW < T >::read (std::istream & s)
  /*!
    Calls to Sequence::operator>> into objects of type ClustalW<T>
    results in a call to this function, which reads the alignment in
    from the stream.
  */
  {
    std::string clustalw;
    char ch;
    std::map < std::string, std::string > seqs;
    //the map "order" is used to make sure that sequences
    //are added to the vector<T*> in the order that they
    //appear in the file--necessary because maps implicity
    //sort their contents...
    std::map < std::string, int >order;
    int nseqs = 0;
    s >> clustalw;
    if (clustalw != "CLUSTAL")
      {
        throw badFormat ("Sequence::ClustalW::read() : input stream does not appear to be in CLUSTALW format");
      }
    else
      {
        ReadThroughLine (s);
      }
    std::string temp, temp2;
    while (!s.eof ())
      {
        s.get (ch);
        bool putback = 0;
        if (ch == '\n' || ch == ' ' || ch == '*')
          {
            s.putback (ch);
            ReadThroughLine (s);
            putback=1;
          }
        else
          {
            if(!putback)
              s.putback(ch);
            s >> temp;
            std::map < std::string, std::string >::iterator iter =
              seqs.find (temp);
            if (iter != seqs.end ())
              {
                s >> temp2 >> std::ws;
                seqs[(*iter).first] += temp2;
              }
            else
              {
                s >> temp2 >> std::ws;
                seqs[temp] = temp2;
                order[temp] = nseqs++;
              }
          }
      }

    typename std::vector<T> _data;
    for (int i = 0; i < nseqs; ++i)
      {
        std::map < std::string, std::string >::iterator iter = seqs.begin (),
	  iter_end = seqs.end ();
        bool found = 0;
        while (iter != iter_end)
          {
            if (order[(*iter).first] == i)
              {
                _data.push_back (T((*iter).first,(*iter).second));
                iter = iter_end;
                found = 1;
              }
            if (!found)
              ++iter;
          }
      }
    this->assign(_data.begin(),_data.end());
    return s;
  }

  template<typename T>
  std::ostream & ClustalW<T>::print(std::ostream &s) const
  {
    typename ClustalW<T>::const_iterator i=this->begin(),
      j=this->end();
    size_t len = i->second.length(),k=0;
    s << "CLUSTAL W" << "\n\n";
    while (k<len)
      {
	size_t offset = (k+60 < len) ? k+60 : k+(len-k);
	for(i=this->begin() ; i < j ; ++i)
	  {
	    s << i->first
	      << '\t';
	    std::copy(i->second.begin()+k,
		      i->second.begin()+offset,
		      std::ostream_iterator<char>(s,""));
	    s << '\n';
	  }
	s << '\n';
	k=offset;
      }
    return s;
  }

  template < typename T >
  std::istream & ClustalW < T >::ReadThroughLine (std::istream & s)
  {
    char ch;
    while (s.get (ch))
      {
        if (ch == '\n')
          return s;
      }
    return s;
  }

}

