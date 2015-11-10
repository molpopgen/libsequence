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
  This is a pure virtual base class for dealing with polymorphism data.  It has no
  real utility in and of itself, other than defining the interface
  to derived classes.  The class publicly inherits from 
  std::pair< std::vector<double>, std::vector<std::string> >,
  representing positions and variable site data, respectively.
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
#include <type_traits>
#include <functional>
#include <memory>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/polySiteVector.hpp>
/*! \example PolyTableIterators.cc */
namespace Sequence
{
  struct PolyTableImpl;
  class PolyTable 
  {
  private:
    std::shared_ptr<PolyTableImpl> impl;
  public:
    //! Data type to store site positions
    using pos_container_t = std::vector<double>;
    //! Data type for storing genotypes
    using geno_container_t = std::vector<std::string>;
    
    //typedefs for container types
    //! \brief non-const reference to std::string
    using reference = geno_container_t::reference;
    //! \brief const reference to std::string
    using const_reference = geno_container_t::const_reference;
    //! \brief The size_type for the haplotype vector
    using size_type = geno_container_t::size_type;
    /*!
      \brief non-const iterator to the haplotypes
    */
    using data_iterator = geno_container_t::iterator;
    /*!
      \brief const iterator to the haplotypes
    */
    using const_data_iterator = geno_container_t::const_iterator;
    /*!
      \brief non-const iterator to the positions
    */
    using pos_iterator =  std::vector<double>::iterator;
    /*!
      \brief const iterator to the positions
    */
    using const_pos_iterator = std::vector<double>::const_iterator;
    /*! \brief Const iterator to segregating sites
      Const iterator to segregating sites. The value_type of this
      iterator is const std::pair<double,std::string>, where the 
      double is the position of the segregating site, and the
      string the list of states at the site.  The first character
      in the string corresponds to the state of the first character
      in the PolyTable (i.e. (*this)[0]), etc.
    */
    using const_site_iterator = Sequence::polySiteVector::const_iterator;

    //functions to return iterators
    data_iterator begin();
    data_iterator end();
    const_data_iterator begin() const;
    const_data_iterator end() const;
    const_data_iterator cbegin() const;
    const_data_iterator cend() const;
    pos_iterator pbegin();
    pos_iterator pend();
    const_pos_iterator pbegin() const;
    const_pos_iterator pend() const;
    const_pos_iterator pcbegin() const;
    const_pos_iterator pcend() const;
    const_site_iterator sbegin() const;
    const_site_iterator send() const;
    const_site_iterator scbegin() const;
    const_site_iterator scend() const;

    //constructor types
    explicit PolyTable();
    explicit PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end);
    explicit PolyTable( std::vector<double> && __positions,
			std::vector<std::string> && __data );
    PolyTable(const PolyTable &) = default;
    PolyTable(PolyTable &) = default;
    PolyTable(PolyTable &&) = default;
    virtual ~ PolyTable (void);

    
    std::vector < double > GetPositions (void) const;
    std::vector < std::string > GetData (void) const;

    //operations on the data (non-const)
    // virtual void ApplyFreqFilter(const unsigned & mincount,
    // 				 const bool & haveOutgroup=false,
    // 				 const unsigned & outgroup = 0);
    // virtual void RemoveMultiHits(const bool & skipOutgroup=false,
    // 				 const unsigned & outgroup=0);
    // virtual void RemoveMissing(const bool & skipOutgroup=false,
    // 			       const unsigned & outgroup=0);
    // virtual void RemoveAmbiguous(const bool & skipOutgroup=false,
    // 				 const unsigned & outgroup=0);
    // virtual void Binary (const bool & haveOutgroup = false,
    // 			 const unsigned & outgroup = 0,
    // 			 const bool & strictInfSites = true);

    //operators and implicit typecasts
    virtual bool operator==(const PolyTable &rhs) const;
    virtual bool operator!=(const PolyTable &rhs) const;
    PolyTable & operator=(PolyTable &&) = default;
    PolyTable & operator=(const PolyTable &) = default;


    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
    const_reference operator[] (const size_type & i) const;

    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
    reference operator[] (const size_type & i);
    
    /*!
      \return true if object contains no data, false otherwise
    */
    bool empty() const;
    
    bool assign(PolyTable::const_site_iterator beg,
		PolyTable::const_site_iterator end);

    bool assign( std::vector<double> && __positions,
		 std::vector<std::string> && __data );

    //! \return Sample size
    size_type size (void) const;

    //! \return the i-th position from the PolyTable::positions.
    double position (const std::vector<double>::size_type & i) const;

    //! \return the number of positions (columns)
    unsigned numsites (void) const;

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
    virtual std::ostream & print(std::ostream &h) const = 0;
  };

  /*!
    \ingroup operators
    Allows objects derived from Sequence::PolyTable
    to be read in from streams
  */
  std::istream & operator>> (std::istream & s, PolyTable & c);

  /*!
    \ingroup operators
    Allows objects derived from Sequence::PolyTable
    to be written out to streams
  */  
  std::ostream & operator<< (std::ostream & o, const PolyTable & c);
}
#endif
