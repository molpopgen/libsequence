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

// Code for the -*- C++ -*- namespace Sequence::PolyTableSlice<T>
#include <Sequence/stateCounter.hpp>
#include <algorithm>

namespace Sequence
{
  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const unsigned & window_size_S,
				     const unsigned & step_len )
    : currentSlice(T()),
      windows( std::vector<range>() ),
      windows_begin(windows.begin()),
      windows_end(windows.end())
				   /*!
				     This constructor calculates sliding windows of a fixed number
				     of segregating sites.
				     \param beg A pointer the first segregating site in the data
				     \param end A pointer to one-past-the-last segregating site in the data
				     \param window_size_S The number of segregating sites in each window
				     \param step_len The number of segregating sites by which to "jump" for each new window
    				   */
  {
    process_windows(beg,end,window_size_S,step_len,1.,false,1.);
  }

  template<typename T>
  PolyTableSlice<T>::PolyTableSlice( const PolyTable::const_site_iterator beg,
				     const PolyTable::const_site_iterator end,
				     const unsigned & window_size,
				     const unsigned & step_len,
				     const double & alignment_length,
				     const double & physical_scale)

    : currentSlice(T()),
      windows( std::vector<range>() ),
      windows_begin(windows.begin()),
      windows_end(windows.end())
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
    process_windows(beg,end,window_size,step_len,
		    alignment_length,true,physical_scale);
  }



  template<typename T>  
  void PolyTableSlice<T>::process_windows( const PolyTable::const_site_iterator beg,
					   const PolyTable::const_site_iterator end,
					   const unsigned & window_size,
					   const unsigned & step_len,
					   const double & alignment_length,
					   const bool & is_physical,
					   const double & physical_scale )
  {
    if (window_size > 0 && step_len > 0 )
      {
	if (end>beg)
	  {
	    unsigned nsites = end-beg;
	    if (is_physical)
	      {
		PolyTable::const_site_iterator b=beg,e=end;
		if (alignment_length>0. && physical_scale > 0.)
		  {
		    double currLen=0.,endOfWindow=0.;
		    unsigned i = 0;
		    currLen = double(i*step_len);
		    while(currLen < alignment_length)
		      {
			endOfWindow = currLen + double(window_size);
			unsigned off1=Sequence::SEQMAXUNSIGNED,
			  off2=Sequence::SEQMAXUNSIGNED;
			size_t j=0;
			while ((b+j) < e)
			  {
			    if ( physical_scale*(b+j)->first > currLen && 
				 physical_scale*(b+j)->first <= endOfWindow )
			      {
				if (off1 == Sequence::SEQMAXUNSIGNED)
				  {
				    off1=j;
				    off2=j+1; //we add 1 because off2 will be used to make an "end" iterator
				  }
				else
				  {
				    ++off2;
				  }
				++j;
			      }
			    else ++j;
			    if ( (b+j)>=e || physical_scale*(b+j)->first > endOfWindow )
			      break;
			  }
			if (off1 != Sequence::SEQMAXUNSIGNED)
			  {
			    windows.push_back( std::make_pair( (b+off1),(b+off2) ) );
			  }
			else
			  {
			    windows.push_back( std::make_pair( end,end ) );
			  }
			++i;
			currLen = double(i*step_len);
		      }
		  }
	      }
	    else
	      {
		unsigned jump=0,k=0;
		for (unsigned i=0 ; i<nsites ; i += jump)
		  {
		    //A check is necessary here--
		    //We are trying to guarantee that there are exactly k polymorphisms
		    //in each window.  However, SNP tables can be constructed/manipulated
		    //such that each column is not neccessarily polymorphic.
		    jump=k=0;
		    while(jump < step_len && i+k<nsites)
		      {
			stateCounter counts = std::for_each( (beg+i+k)->second.begin(),
							     (beg+i+k)->second.end(),
							     stateCounter('-'));
			if (counts.nStates() > 1) //is a polymorphic site
			  {
			    ++jump;
			  }
			++k;
		      }
		    windows.push_back( std::make_pair( (beg+i),(beg+i+k) ) );
// 		    if ( (beg+i+step_len+1) < end )
// 		      {
// 			windows.push_back( std::make_pair( (beg+i),(beg+i+step_len+1) ) );
// 		      }
// 		    else
// 		      {
// 			windows.push_back( std::make_pair( (beg+i),end ) );
// 		      }
		  }
	      }
	    windows_begin = windows.begin();
	    windows_end = windows.end();
	  }
      }
  }

  template<typename T>
  typename PolyTableSlice<T>::const_iterator PolyTableSlice<T>::begin() const
  {
    return windows_begin;
  }

  template<typename T>
  typename PolyTableSlice<T>::const_iterator PolyTableSlice<T>::end() const
  {
    return windows_end;
  }

  template<typename T>  
  T PolyTableSlice<T>::get_slice(const_iterator itr) const
    

  /*!
    \return The window pointed to by the iterator itr.
    \

  */
  {
    if (itr >= windows_end)
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
