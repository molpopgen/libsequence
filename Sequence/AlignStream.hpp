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

#ifndef __ALIGNSTREAM_H__
#define __ALIGNSTREAM_H__

#include <Sequence/Alignment.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <utility>
#include <string>
/*! \file AlignStream.hpp
  @brief declaration of virtual base class for alignment streams (Sequence::AlignStream)
*/

/*!
  \class Sequence::AlignStream Sequence/AlignStream.hpp
  \ingroup alignment
 
  This template defines a virtual base class for 
  input of sequences in aligned block formats, such
  as ClustalW.
 
  This class doesn't really define any operations 
  on alignments.  It is merely an interface for input/output
  and some simple basic manipulations.

  @note The only valid template parameters for this
  class are types in the inheritance hierarchy of 
  Sequence::Seq.  std::pair<std::string,std::string>
  is also supported.

  @short Virtual interface to alignment streams
*/

namespace Sequence
  {
  template < typename T >
  class AlignStream
    {
    private:
      //only compile these templates with types in the class hierarchy of Sequence::Seq
      typedef std::pair<std::string,std::string> baseType;
      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<baseType,T>::value
			    || boost::is_same<baseType,T>::value) );
      /*!
      data is where we store the  objects
      of type T in the alignment.
      */
      std::vector < T >data;
    public:
      AlignStream(const std::vector<T> & _data);
      AlignStream(const AlignStream<T> &a)
      {this->assign(a.begin(),a.end());}
      AlignStream(void){}
      virtual ~ AlignStream (void);
      typedef typename std::vector<T>::size_type size_type;
      typedef typename std::vector<T>::reference reference;
      typedef typename std::vector<T>::const_reference const_reference;
      size_type size(void) const
      /*!
        Returns data.size()
      */
        {
          return data.size();
        }
      reference operator[](const size_type & i) 
      /*!
        Returns the i-th object of type T
        in the vector data
      */
      {
        return data[i];
      }
      const_reference operator[](const size_type & i) const
      /*!
        Returns the i-th object of type T
        in the vector data
      */
      {
        return data[i];
      }
      /*!
	value type is std::vector<T>::iterator 
      */
      typedef typename std::vector<T>::iterator iterator;
      /*!
	value type is std::vector<T>::const_iterator 
      */
      typedef typename std::vector<T>::const_iterator const_iterator;
      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;
      bool IsAlignment (void);
      bool Gapped (void);
      unsigned UnGappedLength (void);
      void RemoveGaps (void);
      void RemoveTerminalGaps (void);
      std::vector < T >Trim ( std::vector < int >sites) 
	;
      std::vector < T >TrimComplement ( std::vector < int >sites) 
	;
      const std::vector< T > Data(void);
      /*!
	To define a non-abstract AlignStream, read must be defined
      */
      virtual std::istream & read (std::istream & s) = 0;
      /*!
	To define a non-abstract AlignStream, print must be defined
      */
      virtual std::ostream & print (std::ostream & s) const = 0;

      /*!
	Assign data to object. Since the value type
	for these iterators evaluates to
	std::vector<T>::const_iterator, any vector<T>
	can be the data source
	\exception Sequence::SeqException is thrown
	if all data elements in the range (beg,end]
	are not of the same length
      */
      void assign(const_iterator beg,const_iterator end)
	;
    };


  template < typename T >
  std::istream & operator>> (std::istream & s, AlignStream < T > &c)
  /*!
    \ingroup operators
    Allows objects derived from Sequence::AlignStream to be read in
    from input streams.  The operator works by calling the virtual
    function Sequence::AlignStream::read
  */
  {
    return c.read (s);
  }

  template < typename T >
  std::ostream & operator<< (std::ostream & s,const  AlignStream < T > &c)
  /*!
    \ingroup operators
    Allows objects derived from Sequence::AlignStream to be printed 
    to output streams.  The operator works by calling the virtual
    function Sequence::AlignStream::print
  */
  {
    return c.print (s);
  }

}
#include <Sequence/bits/AlignStream.tcc>
#endif
