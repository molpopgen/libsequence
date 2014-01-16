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

#include <Sequence/SimData.hpp>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace Sequence
{
  SimData::SimData (const size_t & nsam, const size_t & nsnps):
    PolyTable(nsam,nsnps),totsam(nsam)
    /*!
      The constructor needs to know the sample size simulated.
      This is easily obtainted using Sequence::SimParams.  An example
      of this is found in the example file tajd.cc
    */
  {
  }

  SimData::SimData(double *pos, char **sample, int nsam, int S):
    PolyTable(),totsam(nsam)
  {
    PolyTable::assign(pos,S,sample,nsam);
  }

  std::istream & SimData::read (std::istream & stream) 

  /*!
    A call to istream Sequence::operator>> on a object of type SimData
    results in a call to this function.  NOTE: This is the C++ way to read
    in SimData, but it's the slow way
  */
  {
    char ch;
    while(1)
      {
	if (!(stream >> ch))
	  {
	    stream.setstate(std::ios::failbit);
	    return stream;
	  }
        if (ch == ':')
          break;
      }
    unsigned ss;
    stream >> ss;
    std::vector<double> _positions;
    std::vector<std::string> _data;
    if (ss > 0)
      {
	_positions.resize(ss);
        while(1)
          {
	    stream>>ch;
            if (ch == ':')
              break;
          }
        for (unsigned i = 0; i < ss; ++i)
          {
	    stream >> _positions[i];
          }
	std::string temp;
	while(1)
	  {
	    if (! (stream >> temp))
	      break;
	    else if (temp == "//")
	      {
		break;
	      }
	    else
	      {
		_data.push_back(temp);
	      }
	  }
      }
    else if (ss == 0)
      {
        _positions.resize(0);
	_data.resize(0);
      }
    //assign data into base class
    PolyTable::assign(&_positions[0],ss,&_data[0],_data.size());
    return stream;
  }

  std::ostream & SimData::print(std::ostream &o) const
  /*!
    print the data to \a o
  */
  {
    o << "//\n";
    o << "segsites: " << (*this).numsites() << '\n';
    if (this->numsites()>0)
      {
	o << "positions:";
	for(unsigned i = 0 ; i < (*this).numsites() ; ++i)
	  {
	    o <<" "<< (*this).position(i);
	  }
      }
    o<<'\n';
    for(unsigned i = 0 ; i < (*this).size() ; ++i)
      {
        if (i < (*this).size() - 1)
          o << (*this)[i] << '\n';
        else
          o << (*this)[i];
      }
    return o;
  }

  int SimData::fromfile( FILE * openedfile )
  /*!
    In practice, simulation analysis is I/O intensive.  This method
    provies a routine to read in objects of type Sequence::SimData
    using C-style I/O routines, rather than the C++ operator>>.  This can
    result in huge efficiency gains
    \return an integer that is the the return value from fscanf, so
    that you can check for EOF, etc...
  */
  {
    char ch;
    int rv; //return value from fscanf
    while(1)
      {
        rv = fscanf(openedfile,"%c",&ch);
	if (rv == EOF) return rv;
        if (ch == ':')
          break;
      }
    unsigned ss;
    rv = fscanf(openedfile,"%u",&ss);
    if (rv == EOF) return rv;

    std::vector<double> _positions;
    std::vector<std::string> _data;

    if (ss > 0)
      {
	_positions.resize(ss);
        while(1)
          {
            rv=fscanf(openedfile,"%c",&ch);
	    if (rv == EOF) return rv;
            if (ch == ':')
              break;
          }
        for (unsigned i = 0; i < ss; ++i)
          {
            rv=fscanf(openedfile,"%lf",&_positions[i]);
	    if (rv == EOF) return rv;
          }
	char *seq = new char[ss+2];
	while(1)
	  { 
	    rv=fscanf(openedfile,"%s",seq);
	    if (rv == EOF) 
	      {
		//this is a special case:
		//EOF while reading data
		//means data have been read,
		//so we don't return EOF
		rv=1;
		break;
	      }
	    else if ( strcmp(seq,"//") == 0)  
	      {
		break;
	      }
	    else
	      {
		_data.push_back( std::string(seq) );
	      }
          }
	delete [] seq;
      }
    else if (ss == 0)
      {
        _positions.resize(0);
	_data.resize(0);
      }
    //assign data into base class
    PolyTable::assign(&_positions[0],ss,&_data[0],_data.size());
    return rv;
  }
}
