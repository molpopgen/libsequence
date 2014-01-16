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

#include <Sequence/SimpleSNP.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/Portability/StringStreams.hpp>
#include <cstring>
#include <cctype>
#include <cassert>
#include <iostream>

using std::vector;
using std::string;

namespace Sequence
{
  std::istream & SimpleSNP::read (std::istream & s)  

    /*!
      \exception Sequence::badFormat
    */
  {
    unsigned nsam, nsites, i, j;
    char ch;
    if(!(s >> nsam >> nsites))
      throw badFormat("SimpleSNP.cc: file did not start with s nsam nsites");
    if (Diploid)
      nsam *= 2;
    std::vector<double> _positions(nsites);

    for (i = 0; i < nsites; ++i)
      {
        if (! (s >> _positions[i] >> std::ws))
          throw badFormat("SimpleSNP.cc: error in processing site positions");
      }
    std::string outgroup,temp,temp2;
    std::getline(s,temp);
    istr check(temp),check2(temp);
    unsigned nc = 0;
    while (! check.eof() )
      {
	check >> temp2 >> std::ws;
	++nc;
      }
    bool haveOGlabel = (nc == nsites + 1) ? true : false;
    _names.resize(nsam+1);
    if(haveOGlabel)
      check2 >> _names[0];
    else
      _names[0]="anc";
    outgroup.resize(nsites);
    for (i = 0; i < nsites; ++i)
      {
 	if(!(check2 >> ch))
 	  {
 	    throw badFormat("SimpleSNP.cc: error reading in seg. sites");
 	  }
        ch = char(toupper(ch));
        switch (ch)
          {
          case '?':
            outgroup[i] = 'N';
            break;
          default:
            outgroup[i] = ch;
            break;
          }
      }
    outgroup[i] = '\0';
    unsigned numn = 0;
    //count # of ambiguous bases in outgroup
    for (i = 0; i <  outgroup.length(); ++i)
      {
        if (char(std::toupper(outgroup[i])) == 'N')
          ++numn;
      }
    unsigned haveOG = 0;

    std::vector<std::string> _data;
    if(numn == outgroup.length())
      //if all outgroup characters are ambiguous,
      //then we have no outgroup, so we don't include
      //it in the data vector, and so we don't allocate space
      {
        _data.resize(nsam);
        for (i = 0 ; i < nsam ; ++i)
          _data[i].resize(nsites);
      }
    else
      //otherwise,we have an outgroup and need to allocate space
      {
        _data.resize(nsam+1);
        for (i = 0 ; i < nsam+1 ; ++i)
          _data[i].resize(nsites);

        _data[0] = outgroup;
        haveOG=1;
        haveOutgroup = true;
      }

    for (i = 0 + haveOG; i < nsam + haveOG; ++i)
      {
        string name;	//don't store the name
        if(!(s >> name))
          throw badFormat("SimpleSNP.cc: error processing sequences");
        //_names[i-(haveOG+haveOGlabel)] = name;
	_names[i-haveOG+1] = name;
        char *temp = new char[nsites+1];
        char *temp2 = NULL;
        if (Diploid)
          {
            //_names[i-(haveOG+haveOGlabel)+1] = name;
	    _names[i-haveOG+2] = name;
            temp2 = new char[nsites+1];
          }
        for (j = 0; j < nsites; ++j)
          {
            if(!(s >> ch))
              throw badFormat("SimpleSNP.cc: error processing sequenes");
            ch = char(toupper(ch));
            if (Diploid)
              {
                switch (char(std::toupper(ch)))	//use IUPAC ambiguity symbols
                  {
                  case '?':
                    temp[j] = 'N';
                    temp2[j] = 'N';
                    break;
                  case 'M':
                    temp[j] = 'A';
                    temp2[j] = 'C';
                    break;
                  case 'R':
                    temp[j] = 'A';
                    temp2[j] = 'G';
                    break;
                  case 'W':
                    temp[j] = 'A';
                    temp2[j] = 'T';
                    break;
                  case 'S':
                    temp[j] = 'C';
                    temp2[j] = 'G';
                    break;
                  case 'Y':
                    temp[j] = 'C';
                    temp2[j] = 'T';
                    break;
                  case 'K':
                    temp[j] = 'G';
                    temp2[j] = 'T';
                    break;
                  default:
                    temp[j] = ch;
                    temp2[j] = ch;
                  }
              }
            else if (isoFemale)
              {
                //not implemented yet!!
              }
            else
              {
                switch (ch)
                  {
                  case '?':
                    temp[j] = 'N';
                    break;
                  default:
                    temp[j] = ch;
                    break;
                  }
              }
          }
        temp[j] = '\0';
        if (Diploid)
          {
            temp2[j] = '\0';
            _data[i++] = temp;
            _data[i] =  temp2;
            delete [] temp2;
          }
        else
          _data[i] = temp;

        delete [] temp;
      }
    //assign the data to the base class
    if (_data.size() != nsam+haveOG)
      {
	throw (Sequence::badFormat("SimpleSNP::read() -- number of sequences does not match input value"));
      }     
    PolyTable::assign(&_positions[0],_positions.size(),&_data[0],_data.size());
    RemoveInvariantColumns(this);
    return s;
  }

  std::ostream & SimpleSNP::print(std::ostream &o) const
  {
    o << this->size()-this->haveOutgroup <<'\t' <<this->numsites() << '\n';

    for(unsigned i = 0 ; i < this->numsites()-1 ; ++i)
      {
	o << this->position(i) << '\t';
      }
    o << this->position(this->numsites()-1) << '\n';

    if (haveOutgroup == true)
      {
	if(_names.empty())
	  {
	    o << "anc ";
	  }
	else
	  {
	    o << _names[0];
	  }
        for(unsigned i = 0 ; i < this->numsites() ; ++i)
          {
            o << '\t' << (*this)[0][i];
          }
        o << '\n';
        for(unsigned i = 1 ; i < this->size() ; ++i)
          {
	    if (_names.empty())
	      {
		o << "seq"<<i;
	      }
	    else
	      {
		o << _names[i];
	      }
            for(unsigned j = 0 ; j < this->numsites() ; ++j)
              {
                o << '\t' << (*this)[i][j];
              }
            if (i < this->size()-1)
              o << '\n';
          }
      }
    else
      {
	o << "anc ";
        for(unsigned i = 0 ; i < this->numsites() ; ++i)
          {
            o << '\t' << 'N';
          }
	o << '\n';
        for(unsigned i = 0 ; i < this->size() ; ++i)
          {
	    if (_names.empty())
	      {
		o << "seq"<<i;
	      }
	    else
	      o << _names[i];
            for(unsigned j = 0 ; j < this->numsites() ; ++j)
              {
                o << '\t' << (*this)[i][j];
              }
            if (i < this->size()-1)
              o << '\n';
          }
      }
    return o;
  }

  std::string SimpleSNP::label(unsigned i) const
    /*!
      \return the label the i-th individual in the data
    */
  {
    assert(i<_names.size());
    if (haveOutgroup && i==0)
      {
        return string("outgroup");
      }
    return _names[i-unsigned(haveOutgroup)];
  }

  bool SimpleSNP::outgroup(void) const
    /*!
      returns \c true if there is outgroup information,
      \c false otherwise
    */
  {
    return haveOutgroup;
  }
      
  void SimpleSNP::set_outgroup(const bool & b)
    /*
      Tells the object whether or not the first sequence
      is an outgroup sequence
     */
  {
    haveOutgroup = b;
  }

}
