// Code for the -*- C++ -*- namespace Sequence::phylipData<T>

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

#include <Sequence/phylipData.hpp>
#include <string>
#include <iterator>
#include <utility>
#include <algorithm>
#include <cctype>

namespace Sequence
{
  template < typename T >
  std::istream & phylipData < T >::read (std::istream & s)
  {
    std::string notAllowed = {"()[]:;,"}; //characters not allowed in sequence name
    unsigned nsam,nsites;
    s >> nsam >> nsites >> std::ws;
    std::vector<T> _data(nsam);
    std::string name(10,' '),temp;

    for(unsigned i = 0 ; i < nsam ; ++i)
      {
	s.read( &name[0], 10*sizeof(char) ); //A name is 10 characters
	for( auto c : notAllowed )
	  {
	    if(name.find(c) != std::string::npos)
	      {
		throw std::runtime_error("Sequence::phylipData::read -- invalid character found in sequence name");
	      }
	  }
	std::getline(s,temp);
	temp.erase( std::remove_if(temp.begin(),temp.end(),[](const char & __ch){ return std::isspace(__ch);}), temp.end() );
	_data[i] = T(name,temp);
      }
    s >> std::ws;
    while(!s.eof())
      {
	for ( unsigned i = 0 ; !s.eof() && i < nsam ; ++i )
	  {
	    std::getline(s,temp);
	    s >> std::ws;
	    temp.erase( std::remove_if(temp.begin(),temp.end(),[](const char & __ch){ return std::isspace(__ch);}), temp.end() );
	    _data[i].second+=temp;
	  }
      }
    this->assign(std::move(_data));
    return s;
  }

  template < typename T >
  std::ostream & phylipData < T >::print (std::ostream & s) const
  {
    size_t nsam = this->size();
    s << nsam 
      << '\t'
      << (*this)[0].second.size() << '\n';

    for( auto __t = this->begin() ; __t != this->end() ; ++__t )
      {
	if ( __t->first.length() >= 10 ) 
	  {
	    std::copy( __t->first.begin(),
		       __t->first.begin()+10,
		       std::ostream_iterator<char>(s,""));
	  }
	else  //it is too short
	  {
	    std::copy( __t->first.begin(),
		       __t->first.end(),
		       std::ostream_iterator<char>(s,""));
	    for( decltype(__t->first.size()) i = __t->first.size() ; i < 10 ; ++i )
	      {
		s << ' ';
	      }
	  }
	s << __t->second;
	if( __t < this->end() - 1 ) s << '\n';
      }
    return s;
  }

 
template<typename T>
phylipData<T> & phylipData<T>::operator=( const AlignStream<T> & rhs)
/*!
  An "assignment operator" member function.
  If a phylipData object was constructed with a value of 
  0 for the sequence name length, this function
  will set _namelen to the max sequence name contained in \a rhs.
  If the object was constructed with a value k > 0, the 
  sequence name length will remain unchanged.
*/
{
  this->assign(rhs.begin(),rhs.end());
  return *this;
}

}
