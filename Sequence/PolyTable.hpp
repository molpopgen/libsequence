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
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/PolyTableManip.hpp>

/*! \example PolyTableIterators.cc */
namespace Sequence
{
  class PolyTable : public std::pair< std::vector<double>, std::vector<std::string> >
  {
  private:
    //! A PolyTable publicly inherits from std::pair< std::vector<double>, std::vector<std::string> >
    using PolyTableBase = std::pair< std::vector<double>, std::vector<std::string> >;
    mutable Sequence::polySiteVector pv;
    mutable bool non_const_access;
  public:
    //typedefs for container types
    //! \brief non-const reference to std::string
    using reference = std::vector<std::string>::reference;
    //! \brief const reference to std::string
    using const_reference = std::vector<std::string>::const_reference;
    //! \brief The size_type for the haplotype vector
    using size_type = std::vector<std::string>::size_type;
    /*!
      \brief non-const iterator to the haplotypes
    */
    using data_iterator = std::vector<std::string>::iterator;
    /*!
      \brief const iterator to the haplotypes
    */
    using const_data_iterator = std::vector<std::string>::const_iterator;
    /*!
      \brief non-const iterator to the positions
    */
    using pos_iterator =  std::vector<double>::iterator;
    /*!
      \brief const iterator to the positions
    */
    using const_pos_iterator = std::vector<double>::const_iterator;
    /*! \brief Const iterator to segregating sites
      Const iterator to segregating sites. The value type of this
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
    explicit PolyTable(const size_t & nsam = 0, const size_t nsnps = 0);   
    explicit PolyTable(PolyTable::const_site_iterator beg,
		       PolyTable::const_site_iterator end);
    /*!
      Template constructors simplify compatibility with external data sources.
      \param pbeg Pointer to start of positions
      \param pend Pointer to end of positions in "C++-ese", meaning 1 past the last value.
      \param dbeg Pointer to start of data
      \param dend Pointer to end of data

      \note pbeg/pend must have value types convertible to double.  dbed/dend must have value types convertible to std::string
     */
    template<typename double_type,
	     typename string_type>
    explicit PolyTable( const double_type & pbeg,
			const double_type & pend,
			const string_type & dbeg,
			const string_type & dend ) : PolyTableBase(std::vector<double>(pbeg,pend),
								   std::vector<std::string>(dbeg,dend) ),
						     non_const_access(true)
    {
    }
    /*!
      Constructor for data coming from C data structures.
      \param pbeg Pointer to start of positions
      \param pend Pointer to end of positions in "C++-ese", meaning 1 past the last value.
      \param __data Data matrix
      \param nsam Number of haplotypes in __data

      \note pbeg/pend must have value types convertible to double. 

      Rough guide to usage:
      double * pos = new double[npos];
      char ** dmatrix = new char *[nsam];
      //fill pos and dmatrix with the right stuff...
      //This assignment fails b/c PolyTable is pure virtual class, but you'll get the idea
      PolyTable p(pos,pos+npos,dmatrix,nsam);
     */
    template<typename double_type>
    explicit PolyTable( const double_type & pbeg,
			const double_type & pend,
			const char ** __data,
			const size_t & nsam ) : PolyTableBase(std::vector<double>(pbeg,pend),
							      std::vector<std::string>(nsam)),
						non_const_access(true)
    {
      for( size_t i = 0 ; i < nsam ; ++i )
	{
	  second[i] = std::string( __data[i] );
	}
    }

    explicit PolyTable( std::vector<double> && __positions,
			std::vector<std::string> && __data );
    PolyTable(PolyTable &) = default;
    PolyTable(PolyTable &&) = default;
    virtual ~ PolyTable (void);

    
    std::vector < double > GetPositions (void) const;
    std::vector < std::string > GetData (void) const;

    //operations on the data (non-const)
    virtual void ApplyFreqFilter(const unsigned & mincount,
				 const bool & haveOutgroup=false,
				 const unsigned & outgroup = 0);
    virtual void RemoveMultiHits(const bool & skipOutgroup=false,
				 const unsigned & outgroup=0);
    virtual void RemoveMissing(const bool & skipOutgroup=false,
			       const unsigned & outgroup=0);
    virtual void RemoveAmbiguous(const bool & skipOutgroup=false,
				 const unsigned & outgroup=0);
    virtual void Binary (const bool & haveOutgroup = false,
			 const unsigned & outgroup = 0,
			 const bool & strictInfSites = true);

    //operators and implicit typecasts
    virtual bool operator==(const PolyTable &rhs) const;
    virtual bool operator!=(const PolyTable &rhs) const;
    PolyTable & operator=(PolyTable &&) = default;
    PolyTable & operator=(const PolyTable &) = default;
    operator Sequence::polySiteVector() const;

    //The functions below are inlined data access routines.
    inline const_reference operator[] (const size_type & i) const
    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
    {
      assert(i<second.size());
      return (second[i]);
    }

    inline reference operator[] (const size_type & i)
    /*!
      Return the i-th element of PolyTable::data.
      \note range-checking done by assert()
    */
    {
      assert(i<second.size());
      non_const_access=true;
      return (second[i]);
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

    bool assign( std::vector<double> && __positions,
		 std::vector<std::string> && __data );

//     template<typename iterator>
//     bool rear_insert( const iterator beg,
// 		      const iterator end );

    inline size_type size (void) const
      /*!
        Return how many std::strings are stored
        in PolyTable::data.
      */
    {
      return second.size();
    }

    inline double position (const std::vector<double>::size_type & i) const
      /*!
        Return the i-th position from the PolyTable::positions.
        \note range-checking done by assert()
      */
    {
      assert( i < first.size());
      return first[i];
    }

    inline unsigned numsites (void) const
      /*!
        Return how many positions are stored in
        PolyTable::positions
      */
    {
      return unsigned(first.size());
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
#include <Sequence/bits/PolyTable.tcc>
#endif
