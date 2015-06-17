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
#include <functional>

namespace Sequence
{
  Fasta::Fasta() : Seq() {}

  Fasta::Fasta (const Seq & seq) : Seq(seq) 
    /*! copy constructor */
  {}
 
  Fasta::Fasta( Seq && seq ) : Seq(std::move(seq))
  {
  }

  std::istream & Fasta::read (std::istream & stream) 
  {
    first.clear();
    second.clear();
    std::string temp;
    int ch = stream.peek();
    if( stream.eof() ) { return stream; }
    if (char(ch) != '>')
      {
        throw badFormat("Fasta.cc: error, file not in FASTA format");
      }
    //Read in name
    //stream >> ch >> std::ws;
    ch = stream.get();
    std::getline(stream,first);
    stream >> std::ws;
    second.reserve(1000);
    while( char( ch = stream.peek() ) != '>' && ! stream.eof() )
      {
	std::getline(stream,temp);
	second += temp;
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
