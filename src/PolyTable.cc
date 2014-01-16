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

#include <Sequence/PolyTable.hpp>
#include <Sequence/stateCounter.hpp>
#include <iostream>
#include <cctype>
/*! \defgroup popgen Molecular Population Genetics
 */
/*!
  \defgroup polytables Classes Related to Polymorphism tables
  \ingroup popgen
*/

namespace Sequence
{ 
  PolyTable::PolyTable(const size_t & nsam, const size_t nsnps ) 
  /*!
    Default constructor
    \param nsam when non-zero, the table is allocated to contain a vector of nsam
    strings, each of length nsnps.  Each string is filled with a blank space (the ' ' character).
    \param nsps when non-zero, the table is allocated to contain nsnps position (each with position 0)
  */
    : positions( std::vector<double>(nsnps,0.) ),
      data( std::vector< std::string >( nsam, std::string(nsnps,' ') ) ),
      pv(Sequence::polySiteVector()),
      non_const_access(true)
  {}

  PolyTable::PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end) : 
    pv(Sequence::polySiteVector()),
    non_const_access(true)
  {
    if (beg>=end)
      {
	return;
      }
    assign(beg,end);
  }

  PolyTable::~PolyTable (void)
  {
  }

  bool PolyTable::empty() const
  {
    return data.empty()||positions.empty();
  }

  bool PolyTable::assign(PolyTable::const_site_iterator beg,
			 PolyTable::const_site_iterator end) 
  /*!
    Assignment operation, allowing a range of polymorphic sites
    to be assigned to a polymorphism table.  This exists mainly
    for two purposes. One is the ability to assign tables from 
    "slices" of other tables.  Second is to facilitate the 
    writing of "sliding window" routines.
    \return true if the assignment was successful, false otherwise.
    The only case where false is returned is if the number of individuals
    at each site is not the constan from beg to end.
   */
  {
    non_const_access = true; //set to true in case an exception is thrown
    positions.clear();
    data.clear();
    if(std::distance(beg,end) < 1) return true;
    positions.resize(end-beg);
    pv.resize(end-beg);
    data.resize(beg->second.length());
    size_t nsam = beg->second.length(),i=0,j=0;
    std::string::const_iterator sb,se;
    while((beg+i)<end)
      {
	pv[i]=*(beg+i);
	if ((beg+i)->second.length() != nsam)
	  {
	    //If we toss an exception, let's make sure we leave an empty object.
	    positions.clear();
	    data.clear();
	    return false;
	  }
	positions[i] = (beg+i)->first;
	sb = (beg+i)->second.begin();
	se = (beg+i)->second.end();
	j = 0;
	while( (sb+j) < se )
	  {
	    data[j] += *(sb+j);
	    ++j;
	  }
	++i;
      }
    non_const_access = false;  //everything worked, all private data are assigned, so set to false
    return true;
  }

  bool PolyTable::operator==(const PolyTable &rhs) const
  {
    if (positions.size() != rhs.positions.size()
        || data.size() != rhs.data.size())
      return false;

    for(unsigned i = 0 ; i < positions.size() ; ++i)
      if (positions[i] != rhs.positions[i])
        return false;

    for(unsigned i = 0 ; i < data.size() ; ++i)
      if (data[i] != rhs.data[i])
        return false;

    return true;
  }

  bool PolyTable::operator!=(const PolyTable &rhs) const
  {
    return !(*this==rhs);
  }

  PolyTable::operator Sequence::polySiteVector() const
  /*!
    allow (implicit) typecast of Sequence::PolyTable to Sequence::polySiteVector
  */
  {
    if(non_const_access==true)
      {
	pv = Sequence::rotatePolyTable(this);
	non_const_access=false;
      }
    return pv;
  }

  PolyTable::data_iterator PolyTable::begin()
  /*!
    \return an iterator pointing to the beginning of the std::vector<string> containing
    the data
  */
  {
    non_const_access=true;
    return data.begin();
  }

  PolyTable::data_iterator PolyTable::end()
  /*!
    \return an iterator pointing to the end of the std::vector<string> containing
    the data
  */
  {
    non_const_access=true;
    return data.end();
  }

  PolyTable::const_data_iterator PolyTable::begin() const
  /*!
    \return a const iterator pointing to the beginning of the std::vector<string> containing
    the data
  */
  {
    return data.begin();
  }

  PolyTable::const_data_iterator PolyTable::end() const
  /*!
    \return a const iterator pointing to the end of the std::vector<string> containing
    the data
  */
  {
    return data.end();
  }

  PolyTable::pos_iterator PolyTable::pbegin()
  /*!
    \return an iterator pointing to the beginning of the list of positions
  */
  {
    non_const_access=true;
    return positions.begin();
  }

  PolyTable::pos_iterator PolyTable::pend()
  /*!
    \return an iterator pointing to the end of the list of positions
  */
  {
    non_const_access=true;
    return positions.end();
  }

  PolyTable::const_pos_iterator PolyTable::pbegin() const
  /*!
    \return a const iterator pointing to the beginning of the list of positions
  */
  {
    return positions.begin();
  }

  PolyTable::const_pos_iterator PolyTable::pend() const
  /*!
    \return a const iterator pointing to the beginning of the list of positions
  */
  {

    return positions.end();
  }

  PolyTable::const_site_iterator PolyTable::sbegin() const
  /*!
    \return an object of type Sequence::PolyTable::const_site_iterator
    These iterators allow access to the columns (segregating sites) of
    polymorphism tables
  */
  {
    if(non_const_access == true)
      {
	pv = Sequence::rotatePolyTable(this);
	non_const_access=false;
      }
    return pv.begin();
  }
  
  PolyTable::const_site_iterator PolyTable::send() const
  /*!
    \return an object of type Sequence::PolyTable::const_site_iterator
    These iterators allow access to the columns (segregating sites) of
    polymorphism tables
  */
  {
    if(non_const_access == true)
      {
	pv = Sequence::rotatePolyTable(this);
	non_const_access=false;
      }
    return pv.end();
  }

  void PolyTable::ApplyFreqFilter(unsigned mincount,bool haveOutgroup,
                                  unsigned outgroup)
  /*!
    go through the data and remove all positions where there is a
    variant at count (# of occurences in the sample) < minfreq
    \param mincount minimum count of a variant in the data.  
    Variants that occur < mincount times are thrown out.
    \param haveOutgroup \c true if an outgroup is present in the data,
    \c false otherwise
    \param outgroup the index in the data array containing the outgroup
    (if present)
  */
  {
    std::vector<double> newpos;
    std::vector<std::string> newdata(data.size());

    for(unsigned i = 0 ; i < positions.size() ; ++i)
      {
        stateCounter Counts;
        for (unsigned j = 0 ; j < data.size() ; ++j)
          {
            if(haveOutgroup==false||
	       (haveOutgroup==true && j != outgroup))
              {
                Counts(data[j][i]);
              }
          }
        bool freq = true;
        if (Counts.a > 0 && Counts.a < mincount)
          {
            freq = false;
          }
        else if (freq == true && Counts.g > 0 && Counts.g < mincount)
          {
            freq = false;
          }
        else if (freq == true && Counts.c > 0 && Counts.c < mincount)
          {
            freq = false;
          }
        else if (freq == true && Counts.t > 0 && Counts.t < mincount)
          {
            freq = false;
          }
        else if (freq == true && Counts.zero > 0 && Counts.zero < mincount)
          {
            freq = false;
          }
        else if (freq == true && Counts.one > 0 && Counts.one < mincount)
          {
            freq = false;
          }
        if(freq == true)
          {
            newpos.push_back(positions[i]);
            for(unsigned j = 0 ; j < data.size() ; ++j)
              {
                newdata[j] += data[j][i];
              }
          }
      }
    //Assign takes care of setting non_const_access = true
    assign(&newpos[0],newpos.size(),&newdata[0],newdata.size());
  }

  void PolyTable::RemoveMultiHits(bool skipOutgroup,
                                  unsigned outgroup)
  /*!
    go through the data and remove all the sites with more
    than 2 states segregating.  By default, this routine also removes sites
    where there are 2 states segregating in the ingroup. and the 
    outgroup (if present) has a 3rd state.
    \param skipOutgroup default is \c false.  If \c true, the character
    state of the outgroup is ignored.
    \param outgroup the index of the outgroup in the data vector
  */
  {
    std::vector<double> newpos;
    std::vector<std::string> newdata(data.size());

    for(unsigned i = 0 ; i < positions.size() ; ++i)
      {
        stateCounter Counts;
        for (unsigned j = 0 ; j < data.size() ; ++j)
          {
            if((skipOutgroup==false)||
	       (skipOutgroup==true && j != outgroup))
              {
                Counts(data[j][i]);
              }
          }
        if(Counts.nStates() <= 2)
          {
            newpos.push_back(positions[i]);
            for(unsigned j = 0 ; j < data.size() ; ++j)
              {
                newdata[j] += data[j][i];
              }
          }
      }
    //Assign takes care of setting non_const_access = true
    assign(&newpos[0],newpos.size(),&newdata[0],newdata.size());
  }

  void PolyTable::RemoveMissing(bool skipOutgroup,
                                unsigned outgroup)
  /*!
    go through the data and remove all the sites with 
    missing data (the character N).
    \param skipOutgroup default is \c false.  If \c true, the character
    state of the outgroup is ignored.
    \param outgroup the index of the outgroup in the data vector
  */
  {
    std::vector<double> newpos;
    std::vector<std::string> newdata(data.size());

    for(unsigned i = 0 ; i < positions.size() ; ++i)
      {
        bool hasMissing=false;
        for (unsigned j = 0 ; j < data.size() ; ++j)
          {
            if((skipOutgroup==false)||
	       (skipOutgroup==true && j != outgroup))
              {
                switch(char(std::toupper(data[j][i])))
                  {
                  case 'N':
                    hasMissing=true;
                    break;
                  }
              }
          }
        if(hasMissing==false)
          {
            newpos.push_back(positions[i]);
            for(unsigned j = 0 ; j < data.size() ; ++j)
              {
                newdata[j] += data[j][i];
              }
          }
      }
    //assign takes care of setting non_const_access = true
    assign(&newpos[0],newpos.size(),&newdata[0],newdata.size());
  }

  void PolyTable::RemoveAmbiguous(bool skipOutgroup,
				  unsigned outgroup)
  /*!
    go through the data and remove all the sites with 
    states other than {A,G,C,T,N,-}
    \param skipOutgroup default is \c false.  If \c true, the character
    state of the outgroup is ignored.
    \param outgroup the index of the outgroup in the data vector
  */
  {
    std::vector<double> newpos;
    std::vector<std::string> newdata(data.size());

    for(unsigned i = 0 ; i < positions.size() ; ++i)
      {
	stateCounter c;
        for (unsigned j = 0 ; j < data.size() ; ++j)
          {
            if((skipOutgroup==false)||
	       (skipOutgroup==true && j != outgroup))
              {
		c(data[j][i]);
              }
          }
	if (c.ndna == 0)
	  {
	    newpos.push_back(positions[i]);
            for(unsigned j = 0 ; j < data.size() ; ++j)
              {
                newdata[j] += data[j][i];
              }
	  }
      }
    //assign takes care of setting non_const_access = true
    assign(&newpos[0],newpos.size(),&newdata[0],newdata.size());
  }
  
  void
  PolyTable::Binary (bool haveOutgroup, unsigned outgroup, bool strictInfSites )
  /*!
    Recode the polymorphism table in 0,1 (binary notation)
    \param haveOutgroup use \c true if an outgroup is present, \c false otherwise
    \param outgroup the index of the outgroup in the data vector used to construct the object
    \param strictInfSites if \c true, 

    \note if haveOutgroup== \c true, then 0 means an ancestral state and 1 a derived state 
    in the resulting.  
    /note If haveOutgroup == true, and there are sites with missing data in the outrgroup
    sequence, those sites are removed from the data, since its assumed you actually want
    to know ancestral/derived for every site
  */
  {
    unsigned nsites = this->numsites();

    std::vector<double> newpositions;
    std::vector<unsigned> _pos_indexes;
    std::vector<std::string> newdata((*this).size());

    std::string _outgroup;
    if(haveOutgroup == true)
      {
        _outgroup = (*this)[outgroup];
      }
    else
      {
        _outgroup = (*this)[0];
      }

    for(unsigned j =0 ; j < nsites ; ++j)
      {
        stateCounter Counts;
        for(unsigned i = 0 ; i < (*this).size() ; ++i)
          {
            Counts ((*this)[i][j]);
          }
        if(Counts.nStates()==2)
          {
            newpositions.push_back(this->position(j));
            _pos_indexes.push_back(j);
          }
      }

    for(unsigned i = 0 ; i < (*this).size() ; ++i)
      {
        for(unsigned j = 0 ; j < _pos_indexes.size() ; ++j)
          {
            if(haveOutgroup==false ||
	       (haveOutgroup==true && _outgroup[j] != 'N'))
              //skip sites where the ancestral information
              //is missing AND an outgroup sequence is present
              {
                if ( (*this)[i][_pos_indexes[j]] == 'N')
                  {
                    //missing data (leave untouched...)
                    newdata[i] += 'N';
                  }
                else if( (*this)[i][_pos_indexes[j]]
                         != _outgroup[_pos_indexes[j]] )
                  {
                    //derived (maybe...)
                    newdata[i] += '1';
                  }
                else if ( (*this)[i][_pos_indexes[j]]
                          == _outgroup[_pos_indexes[j]] )
                  {
                    //ancestral
                    newdata[i] += '0';
                  }

              }
          }
      }
    //assign takes care of setting non_const_access = true
    assign(&newpositions[0],newpositions.size(),&newdata[0],newdata.size());
  }

     
  std::vector < double >
  PolyTable::GetPositions (void) const
  /*!
    Returns PolyTable::positions.
  */
  {
    return positions;
  }

  std::vector <std::string > PolyTable::GetData (void) const
  /*!
    Returns PolyTable::data, a vector of std::strings containing polymorphic
    sites. Assuming the vector is returned to a vector<string>
    called data, accessing data[i][j] accesses the j-th site
    of the i-th sequence
  */
  {
    return data;
  }

}
