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
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/SeqExceptions.hpp>

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
    is calculated for non-overlapping 100bp windows.  It is assumed that 1000bp
    were simulated:
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
    Sequence::PolyTableSlice<Sequence::SimData> windows(d.sbegin(),d.send(),100,100,1000,1000);
    PolyTableSlice<Sequence::SimData>::const_iterator itr = windows.begin();
    while(itr < windows.end()) //iterate over windows
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
  */
namespace Sequence 
{
  template<typename T> class PolyTableSlice
  {
  private:
    BOOST_STATIC_ASSERT( (boost::is_base_and_derived<
			  Sequence::PolyTable,T>::value) );
    mutable T currentSlice;
    typedef std::pair<PolyTable::const_site_iterator,
		      PolyTable::const_site_iterator> range;
    //we store the window info as pointers to the range of sites 
    //in each window
    std::vector< range > windows;
    std::vector< range >::const_iterator windows_begin,windows_end;
    void process_windows(  const PolyTable::const_site_iterator beg,
			   const PolyTable::const_site_iterator end,
			   const unsigned & window_size,
			   const unsigned & step_len,
			   const double & alignment_length,
			   const bool & is_physical,
			   const double & physical_scale = 1.);

  public:
    explicit PolyTableSlice(  const PolyTable::const_site_iterator beg,
			      const PolyTable::const_site_iterator end,
			      const unsigned & window_size_S,
			      const unsigned & window_step_len );

    explicit PolyTableSlice(  const PolyTable::const_site_iterator beg,
			      const PolyTable::const_site_iterator end,
			      const unsigned & window_size,
			      const unsigned & step_len,
			      const double & alignment_length,
			      const double & physical_scale = 1.);
    /*!
      const_iterator type to access windows
    */
    typedef std::vector< std::pair<PolyTable::const_site_iterator,
				   PolyTable::const_site_iterator> >::const_iterator const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
    T get_slice(const const_iterator) const ;
    unsigned size() const;
    T operator[](const unsigned &) const ;
  };
}
#include <Sequence/bits/PolyTableSlice.tcc>

#endif
