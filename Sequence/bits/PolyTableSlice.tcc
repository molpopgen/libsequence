// Code for the -*- C++ -*- namespace Sequence::PolyTableSlice<T>

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

#include <Sequence/stateCounter.hpp>
#include <algorithm>
#include <iostream>
namespace Sequence
{
  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const unsigned & window_size_S,
				     const unsigned & step_len )
    : currentSlice(T()),
      windows( std::vector<range>() )
      /*!
	This constructor calculates sliding windows of a fixed number
	of segregating sites.
	\param beg A pointer the first segregating site in the data
	\param end A pointer to one-past-the-last segregating site in the data
	\param window_size_S The number of segregating sites in each window
	\param step_len The number of segregating sites by which to "jump" for each new window
      */
  {
    process_windows_fixed(beg,end,window_size_S,step_len);
  }

  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const double & window_size,
				     const double & step_len,
				     const double & starting_pos)

    : currentSlice(T()),
      windows( std::vector<range>() )
      /*!
	Use this constructor to generate a sliding window accross the sequence itself.
	\param beg A pointer the first segregating site in the data
	\param end A pointer to one-past-the-last segregating site in the data
	\param window_size The size of the sliding window
	\param step_len The distance by which the window jumps
	\param alignment_length The length of the alignment in base pairs.
	\param physical_scale.  For SNP data, set this to 1.  For data with positions labelled
	on the interval [0,1), set this equal to alignment_length.  For example,
	if you simulate data for a 1000bp region using Hudson's program "ms", set this to 1000.
      */
  {
    process_windows(beg,end,window_size,step_len,starting_pos);
  }


  template<typename T>
  void PolyTableSlice<T>::process_windows_fixed(const PolyTable::const_site_iterator beg,
						const PolyTable::const_site_iterator end,
						const unsigned & window_size_S,
						const unsigned & window_step_len )
  {
    std::vector<PolyTable::const_site_iterator> variable_pos;

    //Step 1: record which sites are actually polymorphic,
    //in case PolyTable contains invariant pos'ns
    auto begc = beg;
    while( begc < end )
      {
	stateCounter counts = std::for_each( begc->second.begin(),
					     begc->second.end(),
					     stateCounter('-'));
	if (counts.nStates() > 1) //is a polymorphic site
	  {
	    variable_pos.push_back(begc);
	  }
	begc++;
      }
    variable_pos.push_back(end);
    if(!variable_pos.empty())
      {
	auto vpitr = variable_pos.begin();
	do
	  {
	    windows.push_back(std::make_pair(*vpitr,std::min(*(vpitr+window_size_S),end)));
	    auto lo = std::distance(vpitr,variable_pos.end());
	    if(lo >= window_step_len && lo > window_size_S)
	      vpitr += window_step_len;
	    else{
	      vpitr++;
	      windows.push_back(std::make_pair(*vpitr,std::min(*(vpitr+window_size_S),end)));
	      vpitr=variable_pos.end();
	    }
	  }
	while( vpitr < variable_pos.end() );
      }
  }
  
  template<typename T>  
  void PolyTableSlice<T>::process_windows( const PolyTable::const_site_iterator beg,
					   const PolyTable::const_site_iterator end,
					   const double & window_size,
					   const double & step_len,
					   const double & starting_pos )
  {
    double wbeg = starting_pos;

    //Obtain ptr to first element with position >= wbeg
    auto wbeg_itr = std::lower_bound(beg,end,wbeg,
				     [](const polymorphicSite & __p, const double & __value){
				       return __p.first < __value;
				     });
    while(wbeg_itr != end)
      {
	double wend = wbeg + window_size;
	//ptr to first element with position > wend
	auto wend_itr = std::lower_bound(wbeg_itr,end,wend,
					 [](const polymorphicSite & __p, const double & __value){
					   return __p.first <= __value;
					 });
	windows.push_back( std::make_pair(wbeg_itr,wend_itr) );
	//update window begin data:
	//Update the starting position of next window
	wbeg += step_len;
	//Update pointer to first data lement whose pos'n is >= wbeg
	wbeg_itr = std::lower_bound(wbeg_itr,end,wbeg,
				     [](const polymorphicSite & __p, const double & __value){
				       return __p.first < __value;
				     });
      }
  }

  template<typename T>
  typename PolyTableSlice<T>::const_iterator PolyTableSlice<T>::begin() const
  {
    return windows.begin();
  }

  template<typename T>
  typename PolyTableSlice<T>::const_iterator PolyTableSlice<T>::end() const
  {
    return windows.end();
  }

  template<typename T>  
  T PolyTableSlice<T>::get_slice(const_iterator itr) const
  /*!
    \return The window pointed to by the iterator itr.
    \

  */
  {
    if (itr >= windows.end())
      throw(Sequence::SeqException("PolyTableSlice<T>::get_slice() -- iterator out of range"));
    if(itr->first != itr->second)
      {
	currentSlice.assign(itr->first,itr->second);
	return currentSlice;
      }
    return T();
  }

  template<typename T>
  unsigned PolyTableSlice<T>::size() const
  /*!
    \return The number of windows stored
  */
  {
    return windows.size();
  }

  template<typename T>
  T PolyTableSlice<T>::operator[](const unsigned & i) const 
  /*!
    \return the i-th window
    \throw Sequence::SeqException if subscript i is out of range
  */
  {
    if (i > windows.size())
      throw(Sequence::SeqException("PolyTableSlice::operator[] -- subscript out of range"));
    if(windows[i].first != windows[i].second)
      {
	currentSlice.assign(windows[i].first,windows[i].second);
	return currentSlice;
      }
    return T();
  }
}
