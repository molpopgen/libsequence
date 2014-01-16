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

#ifndef POLYTABLE_H
#define POLYTABLE_H
/*! \file PolyTable.hpp
  @brief Sequence::PolyTable, a virtual base class for polymorphism tables
*/
/*! \class Sequence::PolyTable Sequence/PolyTable.hpp
  \ingroup polytables
  This is a base class for dealing with polymorphism data.  It has no
  real utility in and of itself, other than defining the interface
  to derived classes
  \note Segregating site positions are stored as double rather than int.
  This is because they can be represented as int (for example, in the case of
  Sequence::PolySites). However, it is also reasonable that positions be 
  described as falling along a continuous interval (as in
  the case of Sequence::SimData or Sequence::MS_Interface).

  This class is case-insensitive.  That is, your data can
  include the characters {A,a,G,g,C,c,T,t,N,n,-,0,1}, and still
  work with this library.
  @short The base class for polymorphism tables
*/

#include <string>
#include <vector>
#include <cassert>
#include <iosfwd>
#include <exception>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/PolyTableManip.hpp>

/*! \example PolyTableIterators.cc */
namespace Sequence
{
  class PolyTable
  {
  private:
    /*!
      \c positions is a std::vector of doubles that
      constains the positions of the variable sites,
      from wherever the aligment began. The positions
      are indexed starting from 1, not zero.
    */
    std::vector<double> positions;
    /*!
      \c data is a std::vector of std::strings representing the
      variable sites themselves.  Each member of the 
      std::vector represents one sequence/haplotype.
    */
    std::vector<std::string> data;
    mutable Sequence::polySiteVector pv;
    mutable bool non_const_access;
  public:
    //typedefs for container types
    typedef std::string & reference;
    typedef const std::string & const_reference;
    typedef std::vector<std::string>::size_type size_type;

    /*!
      non-const iterator to the data
    */
    typedef std::vector<std::string>::iterator data_iterator;
    /*!
      const iterator to the data
    */
    typedef std::vector<std::string>::const_iterator const_data_iterator;
    /*!
      non-const iterator to the positions
    */
    typedef std::vector<double>::iterator pos_iterator;
    /*!
      const iterator to the positions
    */
    typedef std::vector<double>::const_iterator const_pos_iterator;
    /*!
      Const iterator to segregating sites. The value type of this
      iterator is const std::pair<double,std::string>, where the 
      double is the position of the segregating site, and the
      string the list of states at the site.  The first character
      in the string corresponds to the state of the first character
      in the PolyTable (i.e. (*this)[0]), etc.
    */
    typedef Sequence::polySiteVector::const_iterator const_site_iterator;

    //functions to return iterators
    data_iterator begin();
    data_iterator end();
    const_data_iterator begin() const;
    const_data_iterator end() const;
    pos_iterator pbegin();
    pos_iterator pend();
    const_pos_iterator pbegin() const;
    const_pos_iterator pend() const;
    const_site_iterator sbegin() const;
    const_site_iterator send() const;

    //constructor types
    explicit PolyTable(const size_t & nsam = 0, const size_t nsnps = 0);   
    explicit PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end);
    virtual ~ PolyTable (void);

    
    std::vector < double > GetPositions (void) const;
    std::vector < std::string > GetData (void) const;

    //operations on the data (non-const)
    virtual void ApplyFreqFilter(unsigned mincount,bool haveOutgroup=false,
				 unsigned outgroup = 0);
    virtual void RemoveMultiHits(bool skipOutgroup=false,
				 unsigned outgroup=0);
    virtual void RemoveMissing(bool skipOutgroup=false,
			       unsigned outgroup=0);
    virtual void RemoveAmbiguous(bool skipOutgroup=false,
				 unsigned outgroup=0);
    virtual void Binary (bool haveOutgroup = false,
			 unsigned outgroup = 0,
       		 bool strictInfSites = true);

    //operators and implicit typecasts
    virtual bool operator==(const PolyTable &rhs) const;
    virtual bool operator!=(const PolyTable &rhs) const;
    operator Sequence::polySiteVector() const;

    //The functions below are inlined data access routines.
    inline const_reference operator[] (const size_type & i) const
    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
    {
      assert(i<data.size());
      return (data[i]);
    }

    inline reference operator[] (const size_type & i)
    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
    {
      assert(i<data.size());
      non_const_access=true;
      return (data[i]);
    }
    /*!
      \return true if object contains no data, false otherwise
    */
    bool empty() const;
    
    bool assign(PolyTable::const_site_iterator beg,
		PolyTable::const_site_iterator end);
    /*! 
      Assign SNP data to the polymorphism table from a vector/array.
      \param _positions an array representing the positions of the SNPs
      \param _num_positions the number of elements in _positions
      \param _data an array containing the characters for each SNP in 
      each individual
      \param _num_individuals the number of elements in _data
      \note If the length of the elements in _data does not
      equal _num_positions, the assignment will fail and you
      will be left with an empty polymorphism table.
      The following piece of code shows how to assign from a std::vector:
      \code
      Sequence::PolySites snpTable;
      std::vector<double> positions;
      std::vector<std::string> data;
      //fill positions and data...
      if ( snpTable.assign(&positions[0],positions.size(),&data[0],data.size()) == true )
      {
      //ok
      }
      else
      {
      //assignment failed for some reason...
      }
      \endcode
    */
    template<typename numeric_type,
	     typename string_type>
    bool assign( const numeric_type * _positions, 
		 const size_t & _num_positions,
		 const string_type * _data,
		 const size_t & _num_individuals );

//     template<typename iterator>
//     bool rear_insert( const iterator beg,
// 		      const iterator end );

    inline size_type size (void) const
      /*!
        Return how many std::strings are stored
        in PolyTable::data.
      */
    {
      return data.size();
    }

    inline double position (const std::vector<double>::size_type & i) const
      /*!
        Return the i-th position from the PolyTable::positions.
        \note range-checking done by assert()
      */
    {
      assert( i < positions.size());
      return positions[i];
    }

    inline unsigned numsites (void) const
      /*!
        Return how many positions are stored in
        PolyTable::positions
      */
    {
      return unsigned(positions.size());
    }

    /*!
      read is a pure virtual function.
      Calls to istream & operator>> (istream & s, PolyTable & c)
      act via this routine, which must be defined in all
      derived classes
    */
    virtual std::istream & read(std::istream &h) = 0;

    /*!
      print is a pure virtual function.
      Calls to ostream & operator<<(ostream & s, PolyTable & c)
      act via this routine, which must be defined in all
      derived classes
    */
    virtual std::ostream & print(std::ostream &h) const =0 ;
  };

  inline std::istream & operator>> (std::istream & s, PolyTable & c)
    /*!
      \ingroup operators
      Allows objects derived from Sequence::PolyTable
      to be read in from streams
    */
  {
    return c.read (s);
  }

  inline std::ostream & operator<< (std::ostream & o, const PolyTable & c)
    /*!
      \ingroup operators
      Allows objects derived from Sequence::PolyTable
      to be written out to streams
    */
  {
    return c.print (o);
  }

}
#include <Sequence/bits/PolyTable.tcc>
#endif
