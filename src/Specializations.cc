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

//Code for the -*- C++ -*- Template Specializations for libsequence
#include <Sequence/Alignment.hpp>
//#include <Sequence/PolyTableSlice.hpp>

/*! \file Specializations.cc
  @brief Definitions of template specializations for library functions
*/

namespace Sequence
{
  namespace Alignment
  {    
    template<> bool Gapped(const std::vector<std::string> &data)
    /*!
      specialization for std::string
    */
    {
      for (int i = 0; unsigned (i) < data.size (); ++i)
        //iterate over sequences
	{
	  if( data[i].find('-') != std::string::npos )
	    return true;
	}

      return false;
    }
    
    template <> bool IsAlignment(const std::vector<std::string> &data)
    /*!
      specialization for std::string
    */
    {
      for (int i = 0; unsigned (i) < data.size (); ++i)
        if (data[i].length () != data[0].length ())
          return 0;
      
      return 1;
    }

    template<>
    bool validForPolyAnalysis( std::vector<std::string>::const_iterator beg,
			       std::vector<std::string>::const_iterator end )
    /*!
      specialization for std::string
    */
    {
      while(beg < end)
	{
	  if (std::find_if(beg->begin(),beg->end(),
			   Sequence::invalidPolyChar())
	      != beg->end())
	    {
	      return false;
	    }
	  ++beg;
	}
      return true;
    }

    template<>
    bool validForPolyAnalysis( std::vector<std::string>::iterator beg,
			       std::vector<std::string>::iterator end )
    {
      while(beg < end)
	{
	  if (std::find_if(beg->begin(),beg->end(),
			   Sequence::invalidPolyChar())
	      != beg->end())
	    {
	      return false;
	    }
	  ++beg;
	}
      return true;
    }

    template <>
    unsigned UnGappedLength(const std::vector <std::string>&data) 

    /*!
      specialization for std::string
    */
    {
      bool site_gapped = 0;
      int len = 0;
      if (!IsAlignment(data))
	return Sequence::SEQMAXUNSIGNED;

      for (size_t j = 0; j < data[0].length (); ++j)

        {
          site_gapped = 0;
          for (size_t i = 0;  i < data.size ();  ++i)
            {
              if (data[i][j] == '-')
                {
                  site_gapped = 1;
                  i = data.size();
                }
            }
          if (!(site_gapped))
            ++len;
        }
      return len;
    }

    template<>
    void RemoveGaps (std::vector <std::string> &data)
    /*!
      a specialization for std::string
    */
    {
      size_t i, j;
      size_t length = data[0].length ();
      std::vector < std::string > ungapped_sequences(data.size());
      bool site_is_gapped;

      for (i = 0; i < length; ++i)
        {	//iterate over sites
          for ( j = 0, site_is_gapped = 0;
                j != data.size();  ++j)
            {
              if (data[j][i] == '-')
                {
                  site_is_gapped = 1;
                  j = data.size();
                }
            }
          if (!(site_is_gapped))
            {
              for ( j = 0 ; j != data.size();  ++j)
                ungapped_sequences[j] += data[j][i];
            }
        }

      //redo the data
      data.assign(ungapped_sequences.begin(),ungapped_sequences.end());
    }

    template<>
    void RemoveTerminalGaps (std::vector <std::string>&data)
    /*!
      a specialization for std::string
    */
    {
      size_t i, j, leftmost, rightmost, numUngapped,offset;
      size_t length = data[0].length();	//how much we have to iterate over
      std::vector < std::string > trimmed_sequences;
      size_t size = data.size();

      leftmost = SEQMAXUNSIGNED;
      rightmost = length + 1;	//offset by one b/c its an array...

      //find the leftmost site where all sites in the alignment are ungapped
      for (i = 0; i < length; ++i)
        {	//iterate over sites
          for (numUngapped = 0, j = 0; j != data.size (); ++j)
            {
              if (data[j][i] != '-')
                ++numUngapped;
            }
          if (numUngapped == size)
            {
              leftmost = i;
              i = length + 1;
            }
        }

      //find the rigthmost site where all sites in the alignment are ungapped
      bool exit_condition = false;
      for (i = length - 1; i < data[0].length() && exit_condition == false; --i)
        {
          for (numUngapped = 0, j = 0; j != data.size (); ++j)
            {
              if (data[j][i] != '-')
                ++numUngapped;
            }
          if (numUngapped == size)
            {
              rightmost = i;
              exit_condition = true;
            }
        }

      //now, fill the array of trimmed sequences
      offset = rightmost - leftmost;
      for (j = 0; j != data.size (); ++j)
        trimmed_sequences.push_back (data[j].substr  (leftmost, rightmost-leftmost+1));

      //now, redo the seq array for the current object
      data.assign(trimmed_sequences.begin(),trimmed_sequences.end());
    }

    template <>
    void RemoveFixedOutgroupInsertions( std::vector<std::string> & data,
					unsigned site,
					const unsigned & ref )
    {
      const size_t nsam = data.size()-1;
      if ( site < data[0].length() )
	{
	  unsigned ngap=0;
	  for(unsigned ind=0;ind<data.size();++ind)
	    {
	      if (ind != ref)
		{
		  ngap += (data[ind][site] == '-') ? 1 : 0;
		}
	    }
	  if(ngap==nsam)
	    {
	      for(unsigned ind=0;ind<data.size();++ind)
		{
		  data[ind].erase(site,1);
		}
	      RemoveFixedOutgroupInsertions(data,site,ref);
	    }
	  RemoveFixedOutgroupInsertions(data,site++,ref);
	}
    }

    template<>
    std::vector <std::string>Trim (const std::vector <std::string >&data,
                                   const std::vector <int> &sites)
      

    /*!
      a specialization for std::string
    */
    {
      size_t i, j, numseqs = data.size (), numIntervals =  sites.size ();
      unsigned start, stop;
      std::vector < std::string >trimmedData(numseqs);
      std::vector < std::string > trimmedTemp(numseqs);
      if (sites.empty ())
        {
          throw SeqException ("Sequence::Alignment::Trim(): empty vector of positions passed to function");
        }
      if (numIntervals % 2 != 0)
        {
          throw SeqException ("Sequence::Alignment::Trim(): odd number of positions passed");
        }

      for (i = 0; i < numIntervals; i += 2)
        {
          start = sites[i];
          stop = sites[i + 1];
          for (j = 0; j < numseqs; ++j)
            {
              trimmedTemp[j] += data[j].substr (start, stop - start + 1);
            }
        }
      trimmedData.assign(trimmedTemp.begin(),trimmedTemp.end());
      return trimmedData;
    }

    template<>
    std::vector <std::string>TrimComplement (const std::vector <std::string>&data,
					     const std::vector < int > &sites) 
      

    /*!
      a specialization for std::string
    */
    {
      std::vector < int >newSites;
      size_t i, j, start, stop, numseqs = data.size (), numIntervals = sites.size (), lastval;

      if (sites.empty ())
        {
          throw SeqException ("Sequence::Alignment::TrimComplement(): empty vector of positions passed to function");
        }
      if (numIntervals % 2 != 0)
        {
          throw SeqException ("Sequence::Alignment::TrimComplement(): odd number of positions passed to function");
        }

      std::vector < std::string >trimmedData(numseqs);
      std::vector < std::string > trimmedTemp(numseqs);

      size_t odd_even;
      if (sites[0] == 0)
        {
          for (i = 1; i < numIntervals; ++i)
            {
              odd_even = i+1;
              if (odd_even%2==0.)
                {
                  newSites.push_back (sites[i] +  1);
                }
              else if (odd_even%2!=0.)
                {
                  newSites.push_back (sites[i] -  1);
                }
            }
        }
      else if (sites[0] > 0)
        {
          newSites.push_back (0);
          for (i = 0; i < numIntervals; ++i)
            {
              odd_even = i+1;
              if (odd_even%2==0.)
                {
                  newSites.push_back (sites[i] +  1);
                }
              else if (odd_even%2!=0.)
                {
                  newSites.push_back (sites[i] -  1);
                }
            }
        }

      lastval = newSites[newSites.size () - 1];
      newSites.pop_back ();
      numIntervals = newSites.size ();
      for (i = 0; i < numIntervals; i += 2)
        {
          start = newSites[i];
          stop = newSites[i + 1];
          for (j = 0; j < numseqs; ++j)
            {
              trimmedTemp[j] +=
                data[j].
                substr (start,  stop - start + 1);
            }
        }
      
      for (j = 0; j < numseqs; ++j)
        {
          trimmedTemp[j] +=data[j].substr (lastval);
        }

      trimmedData.assign(trimmedTemp.begin(),trimmedData.end());
      return trimmedData;
    }
  }
}

