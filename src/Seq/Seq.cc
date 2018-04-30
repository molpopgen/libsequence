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
    Seq::Seq(void) : name{}, seq{}
    /*!
  	creates an empty sequence
      */
    {
    }

    std::string
    Seq::substr(std::string::size_type beg, std::string::size_type len) const
    /*!
      Mimics the std::string member function of the same name.
    */
    {
        return seq.substr(beg, len);
    }

    std::string
    Seq::substr(std::string::size_type beg) const
    /*!
      Mimics the standardstd::string member function of the same name.
    */
    {
        return seq.substr(beg);
    }

    Seq::size_type
    Seq::length(void) const
    /*!
      Return the total length of the sequence in bytes
    */
    {
        return seq.length();
    }

    Seq::size_type
    Seq::size(void) const
    /*!
      Return the total length of the sequence in bytes
    */
    {
        return seq.size();
    }

    Seq::reference Seq::operator[](const size_type &i)
    /*!
      Return the i-th element of the sequence.
      \note range-checking is done by assert()
    */
    {
        assert(i < seq.length());
        return seq[i];
    }

    Seq::const_reference Seq::operator[](const size_type &i) const
    /*!
      Return the i-th element of the sequence.
      \note range-checking is done by assert()
    */
    {
        assert(i < seq.length());
        return seq[i];
    }

    bool
    Seq::operator==(const Seq &rhs) const
    /*!
    \return true if the sequences contain the same data,
    false otherwise.
    \note only the sequences (i.e. this->seq and rhs.seq) are compared
  */
    {
        return (this->seq == rhs.seq);
    }

    bool
    Seq::operator!=(const Seq &rhs) const
    /*!
    \return false if the sequences contain the same data,
    true otherwise.
    \note only the sequences (i.e. this->seq and rhs.seq) are compared
  */
    {
        return (this->seq != rhs.seq);
    }

    Seq::operator std::string() const
    /*!
    allows (implict) cast to std::string
  */
    {
        return seq;
    }

    Seq::size_type
    Seq::UngappedLength(void) const
    /*!
      Return length of sequence, excluding the gap character '-'
    */
    {
        size_type ngap
            = size_type(std::count(seq.begin(), seq.end(), '-'));
        return (seq.length() - ngap);
    }

    bool
    Seq::IsGapped(void) const
    /*!
      Returns 1 if the sequence contaings the gap character '-',
      0 otherwise
    */
    {
        return (seq.find('-') != std::string::npos);
    }

    void
    Seq::Subseq(const unsigned &beg, const unsigned &length)
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
        assert(beg < seq.length() && (beg + length < seq.length()));
        seq.assign(seq.begin() + beg, seq.begin() + beg + length);
    }

    void
    Seq::Complement(void)
    /*!
      Complement the Sequence
      \note This modifies the data in the object
      by changing the std::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
    */
    {
        std::for_each(seq.begin(), seq.end(), ComplementBase());
    }

    void
    Seq::Revcom(void)
    /*!
      Reverse and complement the sequence.  
      \note This function modifies the data in the object
      by changing the std::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
      \return *this
    */
    {
        std::reverse(seq.begin(), seq.end());
        std::for_each(seq.begin(), seq.end(), ComplementBase());
    }

    Seq::iterator
    Seq::begin()
    /*!
      \return an iterator to the beginning of the sequence
    */
    {
        return seq.begin();
    }

    Seq::iterator
    Seq::end()
    /*!
      \return an iterator to the end of the sequence
    */
    {
        return seq.end();
    }

    Seq::const_iterator
    Seq::begin() const
    /*!
      \return a const iterator to the beginning of the sequence
    */
    {
        return seq.begin();
    }

    Seq::const_iterator
    Seq::end() const
    /*!
      \return a const iterator to the end of the sequence
    */
    {
        return seq.end();
    }

    Seq::const_iterator
    Seq::cbegin() const
    /*!
      \return a const iterator to the beginning of the sequence
    */
    {
        return seq.cbegin();
    }

    Seq::const_iterator
    Seq::cend() const
    /*!
      \return a const iterator to the end of the sequence
    */
    {
        return seq.cend();
    }

    const char *
    Seq::c_str(void) const
    /*!
      \return the the C-style string representing the 
      sequence as a cont char *
    */
    {
        return seq.c_str();
    }

    //non-member functions
    std::ostream &
    operator<<(std::ostream &s, const Seq &c)
    {
        return c.print(s);
    }

    std::istream &
    operator>>(std::istream &s, Seq &c)
    {
        return c.read(s);
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
