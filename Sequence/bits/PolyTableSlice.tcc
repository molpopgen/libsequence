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
#include <stdexcept>
#include <cmath>
namespace Sequence
{
  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const unsigned & window_size_S,
				     const unsigned & step_len )
    : windows( std::vector<range>() )
  {
    if(!window_size_S)
      throw std::logic_error("window size cannot be 0");
    if(!step_len)
      throw std::logic_error("step_len cannot be 0");
    if(!std::is_sorted(beg,end,[](const polymorphicSite & a, 
				  const polymorphicSite & b)
		       {
			 return a.first<b.first;
		       }))
      {
	throw std::runtime_error("range (beg,end) must be sorted in increasing order");
      }
    process_windows_fixed(beg,end,window_size_S,step_len);
  }

  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const unsigned nwindows)
    : windows( std::vector<range>() )
  {
    if(!std::is_sorted(beg,end,[](const polymorphicSite & a, 
				  const polymorphicSite & b)
		       {
			 return a.first<b.first;
		       }))
      {
	throw std::runtime_error("range (beg,end) must be sorted in increasing order");
      }
    if( end-beg < nwindows )
      {
	std::cerr << "here\n";
	process_windows_fixed(beg,end,1,1);
      }
    else
      {
	unsigned snp_per_window = unsigned(std::ceil(double(end-beg)/double(nwindows)));
	process_windows_fixed(beg,end,snp_per_window,snp_per_window);
      }
  }

  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const double & window_size,
				     const double & step_len,
				     const double & starting_pos,
				     const double & ending_pos)
    : windows( std::vector<range>() )
  {
    if(window_size<=0.)
      throw std::logic_error("window_size must be > 0");
    if(step_len <= 0.)
      throw std::logic_error("step_len must be > 0");
    if(!std::is_sorted(beg,end,[](const polymorphicSite & a, 
				  const polymorphicSite & b)
		       {
			 return a.first<b.first;
		       }))
      {
	throw std::runtime_error("range (beg,end) must be sorted in increasing order");
      }
    process_windows(beg,end,window_size,step_len,starting_pos,ending_pos);
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
    for(auto begc = beg; begc < end ; ++begc )
      {
	stateCounter counts = std::for_each( begc->second.begin(),
					     begc->second.end(),
					     stateCounter('-'));
	if (counts.nStates() > 1) //is a polymorphic site
	  {
	    variable_pos.push_back(begc);
	  }
      }
    variable_pos.push_back(end);
    if(!variable_pos.empty())
      {
	for( auto vpitr = variable_pos.begin() ; vpitr < variable_pos.end() ; vpitr += window_step_len )
	  {
	    windows.emplace_back( std::make_pair(*vpitr, (window_step_len < std::distance(*vpitr,end)) ? *(vpitr+window_size_S) : end) );
	  }
      }
  }
  
  template<typename T>  
  void PolyTableSlice<T>::process_windows( const PolyTable::const_site_iterator beg,
					   const PolyTable::const_site_iterator end,
					   const double & window_size,
					   const double & step_len,
					   const double & starting_pos,
					   const double & ending_pos)
  {
    double wbeg = starting_pos;

    //Obtain ptr to first element with position >= wbeg
    auto wbeg_itr = std::lower_bound(beg,end,wbeg,
				     [](const polymorphicSite & __p, const double & __value){
				       return __p.first < __value;
				     });
    int nwindows=int((ending_pos-wbeg)/step_len);
    for(int win = 0 ; win < nwindows ; ++win)
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
  typename PolyTableSlice<T>::const_iterator PolyTableSlice<T>::cbegin() const
  {
    return windows.begin();
  }

  template<typename T>
  typename PolyTableSlice<T>::const_iterator PolyTableSlice<T>::cend() const
  {
    return windows.end();
  }

  template<typename T>  
  T PolyTableSlice<T>::get_slice(const_iterator itr) const
  {
    if (itr >= windows.end())
      throw(std::out_of_range("PolyTableSlice<T>::get_slice() -- iterator out of range"));
    if(itr->first != itr->second)
      {
	T rv;
	rv.assign(itr->first,itr->second);
	return rv;
      }
    return T();
  }

  template<typename T>
  typename std::vector<typename PolyTableSlice<T>::range>::size_type PolyTableSlice<T>::size() const
  {
    return windows.size();
  }

  template<typename T>
  T PolyTableSlice<T>::operator[](const unsigned & i) const 
  {
    if (i > windows.size())
      throw(std::out_of_range("PolyTableSlice::operator[] -- subscript out of range"));
    if(windows[i].first != windows[i].second)
      {
	T rv;
	rv.assign(windows[i].first,windows[i].second);
	return rv;
      }
    return T();
  }
}
