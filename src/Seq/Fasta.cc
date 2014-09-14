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

#include <Sequence/Fasta.hpp>
#include <stdexcept>
#include <iostream>

namespace Sequence
{
  /*! \example baseComp.cc */
  Fasta::Fasta (const std::string &name, const std::string &seq):
    Seq(name.c_str(),seq.c_str())
	       /*!
		 constructor for const std::string
	       */
  {}
  Fasta::Fasta (const char *name,const char *seq):Seq(name,seq)
	       /*!
		 constructor for const char *
	       */
  {}
  
  Fasta::Fasta (const Seq & seq) : Seq(seq) 
    /*! copy constructor */
  {}

  std::istream & Fasta::read (std::istream & stream) 

  {
    char ch;
    bool seqflag;

    if (!(stream >> ch))
      {
        stream.setstate (std::ios::badbit);
        return (stream);
      }
    else
      stream.putback (ch);

    std::string temp;

    stream >> ch;

    if (ch != '>')
      {
        throw badFormat("Fasta.cc: error, file not in FASTA format");
      }
    first.clear();
    while (1)
      {
        stream.get (ch);
        if (ch == '\n')
          break;
	first += ch;
      }
    seqflag = 1;
    second.clear();
    second.reserve(1000);
    while (seqflag)
      {
        stream >> ch;
        if (ch == '>')
          {
            stream.putback (ch);
            seqflag = 0;
          }
        else if (stream.eof ())
          seqflag = 0;
        else
          {
	    second += ch;
	    std::getline(stream,temp);
	    second += temp;
          }
      }
    return (stream);
  }

  std::ostream & Fasta::print (std::ostream & stream) const
  {
    stream << '>'
	   << first
	   << '\n'
	   << second;
    return stream;
  }
}
