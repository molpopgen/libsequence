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

#include <cassert>
#include <algorithm>
#include <functional>
#include <Sequence/Seq.hpp>
#include <Sequence/SeqFunctors.hpp>

namespace Sequence
{
  Seq::Seq (void) : SeqBase()
      /*!
  	creates an empty sequence
      */
  {
  }
  
  Seq::Seq (const char *name, const char *seq): SeqBase(name,seq)
  {
  }


  Seq::Seq (const std::string & name, const std::string & seq): SeqBase(name,seq)
  {
  }
  
  Seq::Seq (std::string && name, std::string && seq) : 
    SeqBase(std::move(name),std::move(seq))
  {
  }

  std::string Seq::GetName (void) const
  /*!
    Return the sequence name
  */
  {
    return first;
  }

  std::string Seq::GetSeq (void) const
    /*!
      Return the sequence itself
    */
  {
    return second;
  }

  std::string Seq::substr(std::string::size_type beg,std::string::size_type len) const
    /*!
      Mimics the std::string member function of the same name.
    */
  {
    return second.substr(beg,len);
  }

  std::string Seq::substr(std::string::size_type beg) const
    /*!
      Mimics the standardstd::string member function of the same name.
    */
  {
    return second.substr(beg);
  }

  Seq::size_type Seq::length (void) const
    /*!
      Return the total length of the sequence in bytes
    */
  {
    return second.length();
  }

  Seq::size_type Seq::size (void) const
    /*!
      Return the total length of the sequence in bytes
    */
  {
    return second.size();
  }

  Seq::reference
  Seq::operator[] (const size_type & i)
    /*!
      Return the i-th element of the sequence.
      \note range-checking is done by assert()
    */
  {
    assert(i < second.length ());
    return second[i];
  }

  Seq::const_reference
  Seq::operator[] (const size_type & i) const
    /*!
      Return the i-th element of the sequence.
      \note range-checking is done by assert()
    */
  {
    assert(i < second.length ());
    return second[i];
  }
  
  bool Seq::operator==(const Seq & rhs) const
  /*!
    \return true if the sequences contain the same data,
    false otherwise.
    \note only the sequences (i.e. this->second and rhs.second) are compared
  */
  {
    return (this->second==rhs.second);
  }

  bool Seq::operator!=(const Seq & rhs) const
  /*!
    \return false if the sequences contain the same data,
    true otherwise.
    \note only the sequences (i.e. this->second and rhs.second) are compared
  */
  {
    return (this->second!=rhs.second);
  }

  Seq::operator std::string() const
  /*!
    allows (implict) cast to std::string
  */
  {
    return second;
  }

  Seq::size_type
  Seq::UngappedLength (void) const
    /*!
      Return length of sequence, excluding the gap character '-'
    */
  {
    size_type ngap = size_type(std::count(second.begin(),second.end(),'-'));
    return (second.length () - ngap);
  }

  bool Seq::IsGapped (void) const
    /*!
      Returns 1 if the sequence contaings the gap character '-',
      0 otherwise
    */
  {
    return (second.find('-') != std::string::npos);
  }

  void Seq::Subseq (const unsigned & beg, const unsigned & length)
    /*!
      \param beg the index along the sequence at which the substring begins
      \param length the length of the subseq
      Acts via std::string.substr().  
      Note that this modifies the data in the object
      by changing thestd::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
      \note range-checking done by assert()
    */
  {
    assert ( beg < second.length() && (beg+length < second.length()));
    second.assign(second.begin()+beg,second.begin()+beg+length);
  }

  void Seq::Complement(void)
    /*!
      Complement the Sequence
      \note This modifies the data in the object
      by changing the std::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
    */
  {
    std::for_each(second.begin(),second.end(),ComplementBase());
  }

  void Seq::Revcom (void)
    /*!
      Reverse and complement the sequence.  
      \note This function modifies the data in the object
      by changing the std::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
      \return *this
    */
  {
    std::reverse(second.begin(),second.end());
    std::for_each(second.begin(),second.end(),ComplementBase());
  }

  Seq::iterator Seq::begin()
    /*!
      \return an iterator to the beginning of the sequence
    */
  {
    return second.begin();
  }

  Seq::iterator Seq::end()
    /*!
      \return an iterator to the end of the sequence
    */
  {
    return second.end();
  }

  Seq::const_iterator Seq::begin() const
    /*!
      \return a const iterator to the beginning of the sequence
    */
  {
    return second.begin();
  }

  Seq::const_iterator Seq::end() const
    /*!
      \return a const iterator to the end of the sequence
    */
  {
    return second.end();
  }

  Seq::const_iterator Seq::cbegin() const
    /*!
      \return a const iterator to the beginning of the sequence
    */
  {
    return second.cbegin();
  }

  Seq::const_iterator Seq::cend() const
    /*!
      \return a const iterator to the end of the sequence
    */
  {
    return second.cend();
  }

  const char *Seq::c_str(void) const
    /*!
      \return the the C-style string representing the 
      sequence as a cont char *
    */
  {
    return second.c_str();
  }

  //non-member functions
  std::ostream &
  operator<< (std::ostream & s, const Seq & c)
  {
    return c.print (s);
  }

  std::istream &
  operator>> (std::istream & s, Seq &c)
  {
    return c.read (s);
  }
} //ns Sequence

/*! \defgroup data General I/O
 */
/*!
  \defgroup seqio Classes for sequence I/O
  \ingroup data
*/

/*!
  \namespace Sequence
  The entirety of this library is defined in namespace Sequence.  
  @short The namespace in which this library resides
*/
