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
#include <cctype>
#include <algorithm>
#include <iostream>
#include <iterator>
/*! \defgroup popgen Molecular Population Genetics
 */
/*!
  \defgroup polytables Classes Related to Polymorphism tables
  \ingroup popgen
*/

namespace Sequence
{
  struct PolyTableImpl
  {
    std::vector<double> pos;
    std::vector<std::string> data;
    polySiteVector pv;
    bool non_const_access;
    PolyTableImpl() : pos(std::vector<double>()),
		      data(std::vector<std::string>()),
		      pv(polySiteVector()),
		      non_const_access(true)
    {
    }
   
    PolyTableImpl( std::vector<double> && __positions,
		   std::vector<std::string> && __data ) : pos(std::forward<std::vector<double> >(__positions)),
							  data(std::forward<std::vector<std::string> >(__data)),
							  pv(polySiteVector()),
							  non_const_access(true)
    {
      std::for_each(std::begin(data),std::end(data),[this](const std::string & __s) {
	  if(this->pos.size()!=__s.size())
	    {
	      this->pos.clear();
	      this->data.clear();
	      this->pv.clear();
	      throw Sequence::SeqException("PolyTable: number of positions != length of data element");
	    }
	});
    }

    bool empty() const { return pos.empty() && data.empty(); }
    void clear() { pos.clear();data.clear();pv.clear(); }
    bool assign(PolyTable::const_site_iterator beg,
		PolyTable::const_site_iterator end);
    bool assign( std::vector<double> && __positions,
		 std::vector<std::string> && __data );
  };

  bool PolyTableImpl::assign(PolyTable::const_site_iterator beg,
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
    this->clear();
    if(std::distance(beg,end) < 1) return true;
    pos.resize(std::vector<double>::size_type(end-beg));
    pv.resize(std::vector<double>::size_type(end-beg));
    data.resize(beg->second.length());
    size_t nsam = beg->second.length();
    std::string::const_iterator sb,se;
    typedef PolyTable::const_site_iterator::difference_type DTYPE;
    DTYPE i=0,j=0;
    while((beg+i)<end)
      {
	pv[unsigned(i)]=*(beg+i);
	if ((beg+i)->second.length() != nsam)
	  {
	    //If we toss an exception, let's make sure we leave an empty object.
	    this->clear();
	    return false;
	  }
	pos[unsigned(i)] = (beg+i)->first;
	sb = (beg+i)->second.begin();
	se = (beg+i)->second.end();
	j = 0;
	while( (sb+j) < se )
	  {
	    data[unsigned(j)] += *(sb+j);
	    ++j;
	  }
	++i;
      }
    non_const_access = false;  //everything worked, all private data are assigned, so set to false
    return true;
  }

  bool PolyTableImpl::assign( std::vector<double> && __positions,
			      std::vector<std::string> && __data )
  {
    non_const_access = true;
    this->clear();
    std::swap(pos,__positions);
    std::swap(data,__data);

    if( std::find_if( std::begin(data),std::end(data),
		      [this](const std::string __s) {
			return __s.size() != pos.size();
		      } ) != data.cend() )
      {
	this->clear();
	return false;
      }
    return true;
  }
  
  PolyTable::PolyTable() : impl(std::shared_ptr<PolyTableImpl>(new PolyTableImpl()))
  {
  }
  PolyTable::PolyTable( std::vector<double> && __positions,
			std::vector<std::string> && __data ) : impl(std::shared_ptr<PolyTableImpl>(new PolyTableImpl(std::forward<std::vector<double> >(__positions),
														     std::forward<std::vector<std::string> >(__data))))
  {
  }

  PolyTable::PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end) : impl(std::shared_ptr<PolyTableImpl>(new PolyTableImpl())) 
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
    return impl->empty();
  }

  bool PolyTable::operator==(const PolyTable &rhs) const
  /*!
    \return true if *this == rhs, false otherwise
    \warning case-sensitive
   */
  {
    return (this->impl->pos == rhs.impl->pos) &&
      (this->impl->data == rhs.impl->data);
  }

  bool PolyTable::operator!=(const PolyTable &rhs) const
  {
    return !(*this==rhs);
  }

  bool PolyTable::assign( std::vector<double> && __positions,
			      std::vector<std::string> && __data )
  {
    return impl->assign(std::forward<std::vector<double> >(__positions),
		 std::forward<std::vector<std::string> >(__data));
  }

  bool PolyTable::assign(PolyTable::const_site_iterator beg,
			     PolyTable::const_site_iterator end)
  {
    return impl->assign(beg,end);
  }
  
  PolyTable::data_iterator PolyTable::begin()
  /*!
    \return an iterator pointing to the beginning of the std::vector<string> containing
    the data
  */
  {
    impl->non_const_access=true;
    return impl->data.begin();
  }

  PolyTable::data_iterator PolyTable::end()
  /*!
    \return an iterator pointing to the end of the std::vector<string> containing
    the data
  */
  {
    impl->non_const_access=true;
    return impl->data.end();
  }

  PolyTable::const_data_iterator PolyTable::begin() const
  /*!
    \return a const iterator pointing to the beginning of the std::vector<string> containing
    the data
  */
  {
    return impl->data.begin();
  }

  PolyTable::const_data_iterator PolyTable::cbegin() const
  /*!
    \return a const iterator pointing to the beginning of the std::vector<string> containing
    the data
  */
  {
    return impl->data.cbegin();
  }

  PolyTable::const_data_iterator PolyTable::end() const
  /*!
    \return a const iterator pointing to the end of the std::vector<string> containing
    the data
  */
  {
    return impl->data.end();
  }

  PolyTable::const_data_iterator PolyTable::cend() const
  /*!
    \return a const iterator pointing to the end of the std::vector<string> containing
    the data
  */
  {
    return impl->data.cend();
  }

  PolyTable::pos_iterator PolyTable::pbegin()
  /*!
    \return an iterator pointing to the beginning of the list of positions
  */
  {
    impl->non_const_access=true;
    return impl->pos.begin();
  }

  PolyTable::pos_iterator PolyTable::pend()
  /*!
    \return an iterator pointing to the end of the list of positions
  */
  {
    impl->non_const_access=true;
    return impl->pos.end();
  }

  PolyTable::const_pos_iterator PolyTable::pbegin() const
  /*!
    \return a const iterator pointing to the beginning of the list of positions
  */
  {
    return impl->pos.begin();
  }

  PolyTable::const_pos_iterator PolyTable::pend() const
  /*!
    \return a const iterator pointing to the beginning of the list of positions
  */
  {

    return impl->pos.end();
  }

  PolyTable::const_pos_iterator PolyTable::pcbegin() const
  /*!
    \return a const iterator pointing to the beginning of the list of positions
  */
  {
    return impl->pos.cbegin();
  }

  PolyTable::const_pos_iterator PolyTable::pcend() const
  /*!
    \return a const iterator pointing to the beginning of the list of positions
  */
  {

    return impl->pos.cend();
  }

  PolyTable::const_site_iterator PolyTable::sbegin() const
  /*!
    \return an object of type Sequence::PolyTable::const_site_iterator
    These iterators allow access to the columns (segregating sites) of
    polymorphism tables
  */
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(std::move(make_polySiteVector(*this)));
	impl->non_const_access=false;
      }
    return impl->pv.begin();
  }
  
  PolyTable::const_site_iterator PolyTable::send() const
  /*!
    \return an object of type Sequence::PolyTable::const_site_iterator
    These iterators allow access to the columns (segregating sites) of
    polymorphism tables
  */
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(std::move(make_polySiteVector(*this)));
	impl->non_const_access=false;
      }
    return impl->pv.end();
  }

  PolyTable::const_site_iterator PolyTable::scbegin() const
  /*!
    \return an object of type Sequence::PolyTable::const_site_iterator
    These iterators allow access to the columns (segregating sites) of
    polymorphism tables
  */
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(std::move(make_polySiteVector(*this)));
	impl->non_const_access=false;
      }
    return impl->pv.cbegin();
  }
  
  PolyTable::const_site_iterator PolyTable::scend() const
  /*!
    \return an object of type Sequence::PolyTable::const_site_iterator
    These iterators allow access to the columns (segregating sites) of
    polymorphism tables
  */
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(std::move(make_polySiteVector(*this)));
	impl->non_const_access=false;
      }
    return impl->pv.cend();
  }

  // void PolyTable::ApplyFreqFilter(const unsigned & mincount,
  // 				  const bool & haveOutgroup,
  //                                 const unsigned & outgroup)
  // /*!
  //   go through the data and remove all positions where there is a
  //   variant at count (# of occurences in the sample) < minfreq
  //   \param mincount minimum count of a variant in the data.  
  //   Variants that occur < mincount times are thrown out.
  //   \param haveOutgroup \c true if an outgroup is present in the data,
  //   \c false otherwise
  //   \param outgroup the index in the data array containing the outgroup
  //   (if present)
  //   \note This is only a frequency filter based on minor allele counts.
  // */
  // {
  //   std::vector<double> newpos;
  //   std::vector<std::string> newdata(second.size());

  //   for(unsigned i = 0 ; i < first.size() ; ++i)
  //     {
  //       stateCounter Counts;
  //       for (unsigned j = 0 ; j < second.size() ; ++j)
  //         {
  //           if(haveOutgroup==false||
  // 	       (haveOutgroup==true && j != outgroup))
  //             {
  //               Counts(second[j][i]);
  //             }
  //         }
  //       bool freq = true;
  //       if (Counts.a > 0 && Counts.a < mincount)
  //         {
  //           freq = false;
  //         }
  //       else if (freq == true && Counts.g > 0 && Counts.g < mincount)
  //         {
  //           freq = false;
  //         }
  //       else if (freq == true && Counts.c > 0 && Counts.c < mincount)
  //         {
  //           freq = false;
  //         }
  //       else if (freq == true && Counts.t > 0 && Counts.t < mincount)
  //         {
  //           freq = false;
  //         }
  //       else if (freq == true && Counts.zero > 0 && Counts.zero < mincount)
  //         {
  //           freq = false;
  //         }
  //       else if (freq == true && Counts.one > 0 && Counts.one < mincount)
  //         {
  //           freq = false;
  //         }
  //       if(freq == true)
  //         {
  //           newpos.push_back(first[i]);
  //           for(unsigned j = 0 ; j < second.size() ; ++j)
  //             {
  //               newdata[j] += second[j][i];
  //             }
  //         }
  //     }
  //   //take care of case where new data are empty
  //   if( newpos.empty()) newdata.clear();
  //   //Assign takes care of setting non_const_access = true
  //   assign(std::move(newpos),std::move(newdata));
  // }

  // void PolyTable::RemoveMultiHits(const bool & skipOutgroup,
  //                                 const unsigned & outgroup)
  // /*!
  //   go through the data and remove all the sites with more
  //   than 2 states segregating.  By default, this routine also removes sites
  //   where there are 2 states segregating in the ingroup. and the 
  //   outgroup (if present) has a 3rd state.
  //   \param skipOutgroup default is \c false.  If \c true, the character
  //   state of the outgroup is ignored.
  //   \param outgroup the index of the outgroup in the data vector
  // */
  // {
  //   std::vector<double> newpos;
  //   std::vector<std::string> newdata(second.size());

  //   for(unsigned i = 0 ; i < first.size() ; ++i)
  //     {
  //       stateCounter Counts;
  //       for (unsigned j = 0 ; j < second.size() ; ++j)
  //         {
  //           if((skipOutgroup==false)||
  // 	       (skipOutgroup==true && j != outgroup))
  //             {
  //               Counts(second[j][i]);
  //             }
  //         }
  //       if(Counts.nStates() <= 2)
  //         {
  //           newpos.push_back(first[i]);
  //           for(unsigned j = 0 ; j < second.size() ; ++j)
  //             {
  //               newdata[j] += second[j][i];
  //             }
  //         }
  //     }
  //   //Assign takes care of setting non_const_access = true
  //   assign(std::move(newpos),std::move(newdata));
  // }

  // void PolyTable::RemoveMissing(const bool & skipOutgroup,
  //                               const unsigned & outgroup)
  // /*!
  //   go through the data and remove all the sites with 
  //   missing data (the character N).
  //   \param skipOutgroup default is \c false.  If \c true, the character
  //   state of the outgroup is ignored.
  //   \param outgroup the index of the outgroup in the data vector
  // */
  // {
  //   std::vector<double> newpos;
  //   std::vector<std::string> newdata(second.size());

  //   for(unsigned i = 0 ; i < first.size() ; ++i)
  //     {
  //       bool hasMissing=false;
  //       for (unsigned j = 0 ; j < second.size() ; ++j)
  //         {
  //           if((skipOutgroup==false)||
  // 	       (skipOutgroup==true && j != outgroup))
  //             {
  //               switch(char(std::toupper(second[j][i])))
  //                 {
  //                 case 'N':
  //                   hasMissing=true;
  //                   break;
  //                 }
  //             }
  //         }
  //       if(hasMissing==false)
  //         {
  //           newpos.push_back(first[i]);
  //           for(unsigned j = 0 ; j < second.size() ; ++j)
  //             {
  //               newdata[j] += second[j][i];
  //             }
  //         }
  //     }
  //   if(newpos.empty()) newdata.clear();
  //   //assign takes care of setting non_const_access = true
  //   assign(std::move(newpos),std::move(newdata));
  // }

  // void PolyTable::RemoveAmbiguous(const bool & skipOutgroup,
  // 				  const unsigned & outgroup)
  // /*!
  //   go through the data and remove all the sites with 
  //   states other than {A,G,C,T,N,-,0,1}
  //   \param skipOutgroup default is \c false.  If \c true, the character
  //   state of the outgroup is ignored.
  //   \param outgroup the index of the outgroup in the data vector
  // */
  // {
  //   std::vector<double> newpos;
  //   std::vector<std::string> newdata(second.size());

  //   for(unsigned i = 0 ; i < first.size() ; ++i)
  //     {
  // 	stateCounter c;
  //       for (unsigned j = 0 ; j < second.size() ; ++j)
  //         {
  //           if((skipOutgroup==false)||
  // 	       (skipOutgroup==true && j != outgroup))
  //             {
  // 		c(second[j][i]);
  //             }
  //         }
  // 	if (c.ndna == 0)
  // 	  {
  // 	    newpos.push_back(first[i]);
  //           for(unsigned j = 0 ; j < second.size() ; ++j)
  //             {
  //               newdata[j] += second[j][i];
  //             }
  // 	  }
  //     }
  //   if(newpos.empty()) newdata.clear();
  //   //assign takes care of setting non_const_access = true
  //   assign(std::move(newpos),std::move(newdata));
  // }
  
  // void
  // PolyTable::Binary (const bool & haveOutgroup, 
  // 		     const unsigned & outgroup, 
  // 		     const bool & strictInfSites )
  // /*!
  //   Recode the polymorphism table in 0,1 (binary notation)
  //   \param haveOutgroup use \c true if an outgroup is present, \c false otherwise
  //   \param outgroup the index of the outgroup in the data vector used to construct the object
  //   \param strictInfSites if \c true, 

  //   \note if haveOutgroup== \c true, then 0 means an ancestral state and 1 a derived state 
  //   in the resulting.  
  //   /note If haveOutgroup == true, and there are sites with missing data in the outrgroup
  //   sequence, those sites are removed from the data, since its assumed you actually want
  //   to know ancestral/derived for every site
  // */
  // {
  //   unsigned nsites = this->numsites();

  //   std::vector<double> newpositions;
  //   std::vector<unsigned> _pos_indexes;
  //   std::vector<std::string> newdata((*this).size());

  //   std::string _outgroup;
  //   if(haveOutgroup == true)
  //     {
  //       _outgroup = (*this)[outgroup];
  //     }
  //   else
  //     {
  //       _outgroup = (*this)[0];
  //     }

  //   for(unsigned j =0 ; j < nsites ; ++j)
  //     {
  //       stateCounter Counts;
  //       for(unsigned i = 0 ; i < (*this).size() ; ++i)
  //         {
  //           Counts ((*this)[i][j]);
  //         }
  //       if(Counts.nStates()==2)
  //         {
  //           newpositions.push_back(this->position(j));
  //           _pos_indexes.push_back(j);
  //         }
  //     }

  //   for(unsigned i = 0 ; i < (*this).size() ; ++i)
  //     {
  //       for(unsigned j = 0 ; j < _pos_indexes.size() ; ++j)
  //         {
  //           if(haveOutgroup==false ||
  // 	       (haveOutgroup==true && _outgroup[j] != 'N'))
  //             //skip sites where the ancestral information
  //             //is missing AND an outgroup sequence is present
  //             {
  //               if ( (*this)[i][_pos_indexes[j]] == 'N')
  //                 {
  //                   //missing data (leave untouched...)
  //                   newdata[i] += 'N';
  //                 }
  //               else if( (*this)[i][_pos_indexes[j]]
  //                        != _outgroup[_pos_indexes[j]] )
  //                 {
  //                   //derived (maybe...)
  //                   newdata[i] += '1';
  //                 }
  //               else if ( (*this)[i][_pos_indexes[j]]
  //                         == _outgroup[_pos_indexes[j]] )
  //                 {
  //                   //ancestral
  //                   newdata[i] += '0';
  //                 }

  //             }
  //         }
  //     }
  //   //assign takes care of setting non_const_access = true
  //   assign(std::move(newpositions),std::move(newdata));
  // }

  std::vector < double >
  PolyTable::GetPositions (void) const
  /*!
    Returns PolyTableBase.first
  */
  {
    return impl->pos;
  }

  std::vector <std::string > PolyTable::GetData (void) const
  /*!
    Returns PolyTable::data, a vector of std::strings containing polymorphic
    sites. Assuming the vector is returned to a vector<string>
    called data, accessing second[i][j] accesses the j-th site
    of the i-th sequence
  */
  {
    return impl->data;
  }
  
  PolyTable::const_reference PolyTable::operator[] (const size_type & i) const
    {
      return (impl->data[i]);
    }

    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
  PolyTable::reference PolyTable::operator[] (const size_type & i)

    {
      impl->non_const_access=true;
      return (impl->data[i]);
    }

  PolyTable::size_type PolyTable::size (void) const
  {
      return impl->data.size();
  }

  double PolyTable::position (const std::vector<double>::size_type & i) const
  {
    return impl->pos[i];
  }

  unsigned PolyTable::numsites (void) const
    {
      return unsigned(impl->pos.size());
    }
  
  //non-member functions
  std::istream & operator>> (std::istream & s, PolyTable & c)
  {
    return c.read (s);
  }
  
  std::ostream & operator<< (std::ostream & o, const PolyTable & c)
  {
    return c.print (o);
  }
} //ns Sequence

