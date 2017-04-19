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

#ifndef __POLY_TABLE_SLICE_HPP__
#define __POLY_TABLE_SLICE_HPP__
#include <vector>
#include <utility>
#include <type_traits>
#include <Sequence/SeqConstants.hpp>

/*! \file PolyTableSlice.hpp
  Template class for sliding window analysis
 */

/*!\example slidingWindow.cc */
/*!\example slidingWindow2.cc */

  /*! \class Sequence::PolyTableSlice Sequence/PolyTableSlice.hpp
    This class is a simple container to store "sliding windows" along
    an object in the inheritance hierarchy of Sequence::PolyTable.

    Sliding windows are used in population genetics to look at
    variation in levels of diversity along a region.  This class supports
    two simple ways to make such windows.  The first is the slide a window of 
    some length (in base pairs) along your sequence, recording the SNPs in each
    window.  The number of base pairs that you move the window each time is the
    "step length."  The second type of window is to slide a window of a constant
    number of segregating sites along the SNP table.  In the latter case, the step
    length is the number of segregating sites by which to move the beginning of the 
    window each time. The two different constructors for this class correspond
    to these two different window types.

    These two types of window are useful in different contexts, and it's up to
    the user to decide which one s/he wants.  Please note that all this class
    does is facilitate the generation of the windows.  It does not address 
    any of the statistical headaches that arise from sliding window analyses.
    These issues include multiple test correction, non-independence of overlapping
    windows, variation in selective constraing along a sequence, and variation
    in power from window to window with respect to hypothesis testing.

    The user should be aware that the approach used in Kreitman and Hudson (1991)
    "Inferring the evoltionary histories of the Adh and Adh-dup loci in
    Drosophila melanogaster from patterns of polymorphism and divergence." Genetics
    127: 565 describe a clever variant of the sliding window.  They slide along the 
    physical sequence, but keep the number of synonymous/silent sites constant.
    Their procedure mitigates some of the difficulties mentioned above, but it is 
    not implemented here because it relies on having an annotation for the SNP
    table available.

    The user is also referred to Andolfatto, P., J. D. Wall and M. Kreitman, 1999
    "Unusual haplotype structure at the proximal breakpoint of In(2L)t in a 
    natural population of Drosophila melanogaster." Genetics 153:1397-1399,
    which discusses the multiple testing issue.
    
    The following example reads in data from Hudson's program ms.  Tajima's D
    is calculated for non-overlapping windows of size 0.1:
    \code
    #include<Sequence/SimData.hpp>
    #include<Sequence/SimParams.hpp>
    #include<Sequence/PolyTableSlice.hpp>
    #include<Sequence/PolySIM.hpp>
    #include<iostream>
    #include<cstdio>

    int main(int argc, char **argv)
    {
    Sequence::SimParams p;
    std::cin >> p;
    Sequence::SimData d(p.totsam());
    int i;
    while( (i=d.fromstdin()) && i != EOF) //read simulated data from stdin
    {
    Sequence::PolyTableSlice<Sequence::SimData> windows(d.sbegin(),d.send(),0.1,0.1,0);
    //The object only alows const iterations, and libsequence 1.8.5 changed
    //the API to use the "C++11-ese" syntax of cbegin/cend"
    PolyTableSlice<Sequence::SimData>::const_iterator itr = windows.cbegin();
    while(itr < windows.cend()) //iterate over windows
    {
    //create a data object for the current window
    SimData window = windows.get_slice(itr);
    //calculate and print Tajima's D for the window
    PolySIM analyze(&window); 
    std::cout << analyze.TajimasD() << '\t';
    }
    std::cout << std::endl;
    }
    }
    \endcode
    \ingroup popgenanalysis
    \short A container class for "sliding windows" along a polymorphism table
    \note The two constructors are left ambiguous intentionally!  See their documentation below.
    The ambiguity is so that the progammer is forced to thing about which type of slding window to 
    use.
  */
namespace Sequence 
{
  template<typename T> class PolyTableSlice
  {
  public:
    //! Range of a window = [first,second)
    typedef std::pair<PolyTable::const_site_iterator,
		      PolyTable::const_site_iterator> range;
  private:
    static_assert( std::is_base_of<Sequence::PolyTable,T>::value,
		   "T must be derived from Sequence::PolyTable" );
    //we store the window info as pointers to the range of sites 
    //in each window
    std::vector< range > windows;
    void process_windows(  const PolyTable::const_site_iterator beg,
			   const PolyTable::const_site_iterator end,
			   const double & window_size,
			   const double & step_len,
			   const double & starting_pos,
			   const double & ending_pos);

    void process_windows_fixed(  const PolyTable::const_site_iterator beg,
				 const PolyTable::const_site_iterator end,
				 const unsigned & window_size_S,
				 const unsigned & window_step_len );

  public:
    /*!
      This constructor calculates sliding windows of a fixed number
      of segregating sites.
      \param beg A pointer the first segregating site in the data
      \param end A pointer to one-past-the-last segregating site in the data
      \param window_size_S The number of segregating sites in each window
      \param step_len The number of segregating sites by which to "jump" for each new window
      \note In order to use this constructor, you must make sure that the compiler sees unsigned values,
      otherwise compilation will fail with an ambiguity error:
      \code
      PolyTableSlice<PolySites> windows(data,100u,10u);
      \endcode

      \throw std::logic_error if window_size_S or window_step_len == 0
    */
    explicit PolyTableSlice( const PolyTable::const_site_iterator beg,
			     const PolyTable::const_site_iterator end,
			     const unsigned & window_size_S,
			     const unsigned & window_step_len );

    /*!
      Create a specific number of windows with an equal number of segregating sites per window.
      \param beg A pointer the first segregating site in the data
      \param end A pointer to one-past-the-last segregating site in the data
      \param nwindows The desired number of windows.
      \note The intended use of this fxn is to break an interval up into approximately equal-sized 
      chunks.  When end-beg is small relative to nwindows, you will end up with fewer than nwindows
      "slices".  The primary use scenario envisioned for this type of window is downstream 
      parallelization of computation on large PolyTable objects.
     */
    explicit PolyTableSlice( const PolyTable::const_site_iterator beg,
			     const PolyTable::const_site_iterator end,
			     const unsigned nwindows);
    

    /*!
      Use this constructor to generate a sliding window accross the sequence itself.
      \param beg A pointer the first segregating site in the data
      \param end A pointer to one-past-the-last segregating site in the data
      \param window_size The size of the sliding window (in units of physical distance)
      \param step_len The distance by which the window jumps (in units of physical distance)
      \param starting_pos The starting position for your data.
      \param ending_pos.  The last position for the data.
      \note For most situations involving "ms-like" data, and probably for most genomic data,
      a starting_pos of 0 is appropriate.  However, there are situations where your data may be
      something like a segment from the middle of a genome, and your SNP positions are annotated
      with respect to the reference contig.  In that case, starting_pos is best set to the appropriate
      position along the reference, else you will be returned a lot of empty windows should you start from
      0.  Likewise, for "normal" ms runs, ending_pos should be 1.0.  For other scenario, such as "real"
      data, you'll have to set starting_pos and ending_pos to the appropriate values.
      In order to use this constructor, you must make sure that the compiler sees doubles,
      otherwise compilation will fail with an ambiguity error:
      \code
      PolyTableSlice<SimData> windows(data,0.1,0.01);
      \endcode

      \throw std::logic_error if window_size or step_len <= 0.
    */
    explicit PolyTableSlice(  const PolyTable::const_site_iterator beg,
			      const PolyTable::const_site_iterator end,
			      const double & window_size,
			      const double & step_len,
			      const double & starting_pos = 0.,
			      const double & ending_pos = 1.0);
    /*!
      const_iterator type to access windows
    */
    typedef std::vector<range>::const_iterator const_iterator;
    /*!
      \return Const iterator to begin
    */
    const_iterator cbegin() const;
    /*!
      \return Const iterator to end
    */
    const_iterator cend() const;
    /*!
      \param itr An iterator from the current object
      \return The window pointed to by the iterator itr.
      \throw std::out_of_range if itr is out of range
    */
    T get_slice(const const_iterator) const ;
    /*!
      \return The number of windows stored
    */
    std::vector<range>::size_type size() const;
    /*!
      \param i The window to return, 0 <= i < object.size()
      \return the i-th window
      \throw std::out_of_range if i is out of range
    */
    T operator[](const unsigned &) const ;
  };
}
#include <Sequence/bits/PolyTableSlice.tcc>

#endif
