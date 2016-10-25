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
#include <stdexcept>

/*! \defgroup popgen Molecular Population Genetics
 */
/*!
  \defgroup polytables Classes Related to Polymorphism tables
  \ingroup popgen
*/

namespace Sequence
{
  struct PolyTable::PolyTableImpl
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
	      throw std::runtime_error("PolyTable: number of positions != length of data element");
	    }
	});
    }

    
    bool empty() const { return pos.empty() && data.empty(); }
    void clear() { pos.clear();data.clear();pv.clear(); }
    bool assign(PolyTable::const_site_iterator beg,
		PolyTable::const_site_iterator end);
    template<typename postype,typename datatype>
    bool assign_data( postype && __pos,
		      datatype && __data )
    {
      non_const_access=true;
      pos=std::forward<postype>(__pos);
      data=std::forward<datatype>(__data);
      if( std::find_if( std::begin(data),std::end(data),
			[this](const std::string & __s) {
			  return __s.size() != pos.size();
			} ) != data.cend() )
	{
	  this->clear();
	  return false;
	}
      return true;
    }

  };

  bool PolyTable::PolyTableImpl::assign(PolyTable::const_site_iterator beg,
					PolyTable::const_site_iterator end)
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
  /*
  bool PolyTable::PolyTableImpl::assign( std::vector<double> && __positions,
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
  */
  PolyTable::PolyTable() : impl(std::unique_ptr<PolyTableImpl>(new PolyTableImpl()))
  {
  }
  
  PolyTable::PolyTable( std::vector<double> __positions,
			std::vector<std::string> __data ) : impl(std::unique_ptr<PolyTableImpl>(new PolyTableImpl(std::move(__positions),
													       std::move(__data))))
  {
  }

  PolyTable::PolyTable(PolyTable && rhs) : impl(std::unique_ptr<PolyTableImpl>(new PolyTableImpl()))
  {
    std::swap(this->impl,rhs.impl);
  }

  PolyTable::PolyTable(const PolyTable & rhs) : impl(std::unique_ptr<PolyTableImpl>(new PolyTableImpl(*rhs.impl.get())))
  {
  }
  
  PolyTable & PolyTable::operator=(PolyTable && rhs)
  {
    auto x = std::unique_ptr<PolyTableImpl>(new PolyTableImpl());
    std::swap(this->impl,rhs.impl);
    std::swap(x,rhs.impl);
    return *this;
  }

  PolyTable & PolyTable::operator=(const PolyTable & rhs)
  {
    this->impl->pos = rhs.impl->pos;
    this->impl->data = rhs.impl->data;
    return *this;
  }
  
  PolyTable::PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end) : impl(std::unique_ptr<PolyTableImpl>(new PolyTableImpl())) 
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
  {
    return (this->impl->pos == rhs.impl->pos) &&
      (this->impl->data == rhs.impl->data);
  }

  bool PolyTable::operator!=(const PolyTable &rhs) const
  {
    return !(*this==rhs);
  }

  bool PolyTable::assign( const std::vector<double> & __positions,
			  const std::vector<std::string> & __data )
  {
    return impl->assign_data(__positions,__data);
  }

  bool PolyTable::assign( std::vector<double> && __positions,
			  std::vector<std::string> && __data )
  {
    return impl->assign_data(std::move(__positions),
			     std::move(__data));
  }

  bool PolyTable::assign(PolyTable::const_site_iterator beg,
			 PolyTable::const_site_iterator end)
  {
    return impl->assign(beg,end);
  }

  void PolyTable::swap( PolyTable & pt)
  {
    std::swap(this->impl,pt.impl);
  }
  
  PolyTable::data_iterator PolyTable::begin()
  {
    impl->non_const_access=true;
    return impl->data.begin();
  }

  PolyTable::data_iterator PolyTable::end()
  {
    impl->non_const_access=true;
    return impl->data.end();
  }

  PolyTable::const_data_iterator PolyTable::begin() const
  {
    return impl->data.begin();
  }

  PolyTable::const_data_iterator PolyTable::cbegin() const
  {
    return impl->data.cbegin();
  }

  PolyTable::const_data_iterator PolyTable::end() const
  {
    return impl->data.end();
  }

  PolyTable::const_data_iterator PolyTable::cend() const
  {
    return impl->data.cend();
  }

  PolyTable::pos_iterator PolyTable::pbegin()
  {
    impl->non_const_access=true;
    return impl->pos.begin();
  }

  PolyTable::pos_iterator PolyTable::pend()
  {
    impl->non_const_access=true;
    return impl->pos.end();
  }

  PolyTable::const_pos_iterator PolyTable::pbegin() const
  {
    return impl->pos.begin();
  }

  PolyTable::const_pos_iterator PolyTable::pend() const
  {
    return impl->pos.end();
  }

  PolyTable::const_pos_iterator PolyTable::pcbegin() const
  {
    return impl->pos.cbegin();
  }

  PolyTable::const_pos_iterator PolyTable::pcend() const
  {

    return impl->pos.cend();
  }

  PolyTable::const_site_iterator PolyTable::sbegin() const
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(make_polySiteVector(*this));
	impl->non_const_access=false;
      }
    return impl->pv.begin();
  }
  
  PolyTable::const_site_iterator PolyTable::send() const
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(make_polySiteVector(*this));
	impl->non_const_access=false;
      }
    return impl->pv.end();
  }

  PolyTable::const_site_iterator PolyTable::scbegin() const
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(make_polySiteVector(*this));
	impl->non_const_access=false;
      }
    return impl->pv.cbegin();
  }
  
  PolyTable::const_site_iterator PolyTable::scend() const
  {
    if(impl->non_const_access == true)
      {
	impl->pv = polySiteVector(make_polySiteVector(*this));
	impl->non_const_access=false;
      }
    return impl->pv.cend();
  }



  std::vector < double >
  PolyTable::GetPositions (void) const
  {
    return impl->pos;
  }

  std::vector <std::string > PolyTable::GetData (void) const
  {
    return impl->data;
  }
  
  PolyTable::const_reference PolyTable::operator[] (const size_type & i) const
    {
      return (impl->data[i]);
    }

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

