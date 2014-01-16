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

/*! \file phylipData.hpp
  \short Sequence::phylipData -- read in phylip alignments
*/
#ifndef __PHYLIPDATA_HPP__
#define __PHYLIPDATA_HPP__

#include <Sequence/AlignStream.hpp>
#include <utility>
#include <string>
namespace Sequence
{
  template < typename T >
  class phylipData: public AlignStream < T >
		   /*!
		     Input of phylip-format alignments
		    */
  {
  private:
    mutable size_t _namelen;
  public:
    phylipData (const unsigned namelen = 0)
      :AlignStream<T>(),
       _namelen(namelen)
    {}
    phylipData(const std::vector<T> & _data,
	       const unsigned namelen=0)
      :AlignStream<T>(_data),
       _namelen(namelen)
    {
    }
    phylipData (const AlignStream<T> &a) :
      AlignStream<T>(a),_namelen(0)
    {
    }
    ~phylipData(void)
    {}
    std::istream & read (std::istream & s);
    std::ostream & print (std::ostream & s) const;
    phylipData<T> & operator=( const AlignStream<T> & rhs);
  };
}
#include <Sequence/bits/phylipData.tcc>
#endif
