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
#include <Sequence/polySiteVector.hpp>
/*! \example PolyTableIterators.cc */
namespace Sequence
{
  class PolyTable 
  {
  private:
    struct PolyTableImpl;
    std::unique_ptr<PolyTableImpl> impl;
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
    using row_t = geno_container_t::value_type;
    using column_t = const_site_iterator::value_type;
    //functions to return iterators
    //! \return an iterator pointing to the first "haplotype"
    data_iterator begin();
    //! \return an iterator pointing the end of the "haplotypes"
    data_iterator end();
    //! \return a const iterator pointing to the first "haplotype"
    const_data_iterator begin() const;
    //! \return a const iterator pointing the end of the "haplotypes"
    const_data_iterator end() const;
    //! \return a const iterator pointing to the first "haplotype"
    const_data_iterator cbegin() const;
    //! \return a const iterator pointing the end of the "haplotypes"
    const_data_iterator cend() const;
    //! \return iterator to first position
    pos_iterator pbegin();
    //! \return iterator to end of positions
    pos_iterator pend();
    //! \return const iterator to first position
    const_pos_iterator pbegin() const;
    //! \return const iterator to end of positions
    const_pos_iterator pend() const;
    //! \return const iterator to first position
    const_pos_iterator pcbegin() const;
    //! \return const iterator to end of positions
    const_pos_iterator pcend() const;
    //! \return const iterator to first column (position, variants)
    const_site_iterator sbegin() const;
    //! \return const iterator to end of columns
    const_site_iterator send() const;
    //! \return const iterator to first column (position, variants)
    const_site_iterator scbegin() const;
    //! \return const iterator to first column (position, variants)
    const_site_iterator scend() const;

    //constructor types
    explicit PolyTable();
    explicit PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end);
    explicit PolyTable( std::vector<double> __positions,
			std::vector<std::string> __data );
    PolyTable(PolyTable &&);
    PolyTable(const PolyTable &);
    virtual ~ PolyTable (void);

    //! Convenience function to return site positions
    std::vector < double > GetPositions (void) const;
    //! Conventience function to return data.  Each string is a "haplotype".
    std::vector < std::string > GetData (void) const;
    
    //operators and implicit typecasts

    //! Comparison operator.  Case-sensitive
    virtual bool operator==(const PolyTable &rhs) const;
    //! Not-equal operator. Case-sensitive
    virtual bool operator!=(const PolyTable &rhs) const;
    //! Move assignment
    PolyTable & operator=(PolyTable &&);
    //! Copy assignment
    PolyTable & operator=(const PolyTable &);
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
    
    //! \return true if object contains no data, false otherwise
    bool empty() const;

    /*!
      Assignment operation, allowing a range of polymorphic sites
      to be assigned to a polymorphism table.  This exists mainly
      for two purposes. One is the ability to assign tables from 
      "slices" of other tables.  Second is to facilitate the 
      writing of "sliding window" routines.
      \return true if the assignment was successful, false otherwise.
      The only case where false is returned is if the number of individuals
      at each site is not the constant from beg to end.
    */
    bool assign(PolyTable::const_site_iterator beg,
		PolyTable::const_site_iterator end);

    //Swap data with another PolyTable
    void swap( PolyTable & );
    /*!
      Assign data to object
      \return true if successful
    */
    bool assign( const std::vector<double> & __positions,
		 const std::vector<std::string> & __data );
    /*!
      Move data to object
      \return true if successful
    */
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
