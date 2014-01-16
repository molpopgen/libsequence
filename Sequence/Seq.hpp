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

#ifndef SEQ_H
#define SEQ_H
/*! \file Seq.hpp
  @brief class Sequence::Seq, an abstract base class for Sequences
*/

/*! \class Sequence::Seq Sequence/Seq.hpp
  \ingroup seqio
  Abstract interface to sequence objects.  A sequence consists of a name and a sequence,
  both of which are stored as C++ std::strings. Sequence::Seq inherits publicly from 
  std::pair<std::string,std::string> to store the sequence name and data.  The public
  member "first" is the name, and "second" is the data.  These can be used to directly
  access the name and/or data.  Most of the rest of the member functions of this class
  are to make the behavior of Sequence::Seq more "std::string"-like.
  @short Abstract interface to sequence objects
*/

#include <iosfwd>
#include <string>
#include <exception>
#include <utility>
#include <Sequence/SeqExceptions.hpp>

namespace Sequence
  {
    class Seq : public std::pair<std::string,std::string>
    {
      typedef std::pair<std::string,std::string> SeqBase;
    public:
      Seq (void);
      Seq (const char *name, const char *seq);
      Seq (const Seq & seq);
      virtual ~ Seq (){}
      std::string GetName (void) const;
      std::string GetSeq (void) const;
      std::string substr(unsigned beg, unsigned len) const;
      std::string substr(unsigned beg) const;
      /*!
	Iterator to sequence elements.  Iterators 
	access the data, not the seuquence name.
	Value type de-references to char
      */
      typedef std::string::iterator iterator;
      /*!
	Const iterator to sequence elements.  Iterators 
	access the data, not the seuquence name.
	Value type de-references to char
      */
      typedef std::string::const_iterator const_iterator;
      typedef std::string::reference reference;
      typedef std::string::const_reference const_reference;
      typedef std::string::size_type size_type;
      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;
      void Revcom (void);
      void Subseq (unsigned,unsigned);
      void Complement(void);
      size_type length (void) const;
      size_type UngappedLength (void) const;
      bool IsGapped (void) const;
      reference operator[] (const size_type & i);
      const_reference operator[] (const size_type & i) const;
      bool operator==(const Seq & rhs) const;
      bool operator!=(const Seq & rhs) const;
      operator std::string() const;
      const char *c_str(void) const;
      /*!
	read an object of type Sequence::Seq from an istream
      */
      virtual std::istream & read (std::istream & s) = 0;
      /*!
	read an object of type Sequence::Seq from an istream
      */
      virtual std::ostream & print (std::ostream & s) const = 0;
    };

  inline std::ostream &
  operator<< (std::ostream & s, const Seq & c)
  /*!
    \ingroup operators
    Allows objects derived from Sequence::Seq
    to be written to output streams.  This operator
    acts by a call to the virtual funtion Sequence::Seq::print
  */
  {
    return c.print (s);
  }

  inline std::istream &
  operator>> (std::istream & s, Seq &c)
  /*!
    \ingroup operators
    Allows objects derived from Sequence::Seq
    to be read from output streams.  This operator
    acts by a call to the virtual funtion Sequence::Seq::read
  */
  {
    return c.read (s);
  }
}
#endif
