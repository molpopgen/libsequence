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

// Code for the -*- C++ -*- namespace Sequence::Alignment
/*! \file Alignment.tcc
  The static assertion that is used throughout this file
  ensures that the template type parameter T is either
  std::pair<std::string,std::string> or derived from that type,
  which includes Sequence::Seq (see Sequence/Seq.hpp).
  @brief implementation of routines declared in Alignment.hpp
*/
#include <Sequence/Alignment.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/SeqProperties.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <iterator>
#include <algorithm>
namespace Sequence
{
  namespace Alignment
  {
    template < typename T >
    void GetData (std::vector < T >&seqarray,
                  const char *infilename) 
    /*!
      Read objects of type T and put them into the vector seqarray.
      Note that seqarray is not const, so that's where the data go.
      infilename refers to a file from which to read the data. Any
      exceptions that can be thrown when reading the data are not 
      caught here. 
      \n
      Example:\n
      \code
      #include <iostream>
      #include <Sequence/Fasta.hpp>
      #include <Sequence/Alignment.hpp>
      //using namespace std;
      using namespace Sequence;
      int main
      {
      vector<Fasta> data;
      const char *infile = foo;
      try {
      Alignment::GetData(&data,foo);
      } catch (SeqException &e) {
      cerr << "whoa--exception caught"<<endl;
      e.print(cerr);
      cerr << endl;
      }
      }
      \endcode
      \param seqarray a  vector<T> that you want filled
      \param infilename name of file from which fo fill seqarray
      \note if \a infilename is NULL, the function returns, having done nothing
    */
    {
      if(!infilename || infilename == NULL)
	return;

      std::ifstream infile(infilename);
      if (infile)
	{
	  T temp;
	  while (!(infile.eof ()))
	    {
	      infile >> temp;
	      seqarray.push_back(temp);
	    }
	}
    }

    template < typename T >
    std::istream & GetData (std::vector <
			    T >&seqarray,
			    std::istream & input_stream) 
    /*!
      Read objects of type T and put them into the vector seqarray.
      Note that seqarray is not const, so that's where the data go.
      This function is similar to GetData(vector<T> &seqarray,
      const char *infilename), except that it is passed a reference
      to an open input stream, such as a file stream, cin, etc.
      \param seqarray an empty vector<T> that you want filled
      \param input_stream name of istream from which fo fill seqarray
    */
    {
      T temp;
      if (input_stream)
	{
	  while (!(input_stream.eof ()))
	    {
	      input_stream >> temp;
	      seqarray.push_back(temp);
	    }
	}
      return input_stream;
    }

    template < typename T >
    std::istream & ReadNObjects ( std::vector < T > &seqarray,
				  unsigned n,
				  std::istream & input_stream)
    /*!
      Read a fixed number n of objects of type T and put 
      them into the vector seqarray.
      Note that seqarray is not const, so that's where the data go.
      This function is similar to GetData(vector<T> &seqarray,
      const char *infilename), except that it is passed a reference
      to an open input stream, such as a file stream, cin, etc.
      \param seqarray an empty vector<T> that you want filled
      \param n number of objects of type T to read
      \param input_stream name of istream from which fo fill seqarray
      \note the stream is not closed after each read of n records!!
    */
    {
      for(unsigned i = 0 ; i < n ; ++i)
        {
          if (!input_stream || input_stream.eof())
            return input_stream;
          T temp;
          input_stream >> temp;
          seqarray.push_back(temp);
        }
      return input_stream;
    }

    template < typename T >
    void EmptyVector (std::vector< T * > &seqarray)
    /*!
      Free all the memory in seqarray by deleting every objet,
      and resize() seqarray to 0
      \param seqarray the vector<T*> you want emptied
    */
    {
      for(unsigned i = 0 ; i < seqarray.size() ; ++i)
        delete seqarray[i];
      seqarray.resize(0);
    }

    template < typename T >
    bool Gapped (const std::vector < T >&data)
    /*!
      \param data vector<T> containing sequence data
      \return true if the vector contains a gap character ('-')
      , false otherwise.
    */
    {
      BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			   || boost::is_same<std::pair<std::string,std::string>,T>::value));
      for (int i = 0; unsigned (i) < data.size (); ++i)
        //iterate over sequences
	{
	  if( data[i].second.find('-') != std::string::npos )
	    return true;
	}
      return false;
    }

    template < typename T >
    bool IsAlignment (const std::vector < T  >&data)
    /*!
      A vector of sequences/strings is only an alignment if all
      strings are the same length.
      \param data vector<T> to check
    */
    {
       BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			    || boost::is_same<std::pair<std::string,std::string>,T>::value));
      for (int i = 0; unsigned (i) < data.size (); ++i)
        if (data[i].second.length () != data[0].second.length ())
          return 0;

      return 1;
    }

    template<typename Iterator>
    bool validForPolyAnalysis( Iterator beg,
			       Iterator end )
    /*!
      \return true if each element in the range [beg,end) only contains
      characters in the set {A,G,C,T,N,-}, false otherwise
    */
    {
      //ensure that Iterator is an iterator type for a container
      //containing things in the inheritance hierarchy of Sequence::Seq
      //The assertion assumes that std::iterator_traits is specialized for Iterator
      BOOST_STATIC_ASSERT( (boost::is_base_and_derived< 
			    std::pair<std::string,std::string>,
			    typename std::iterator_traits<Iterator>::value_type
			    >::value || 
			    boost::is_same< 
			    std::pair<std::string,std::string>,
			    typename std::iterator_traits<Iterator>::value_type
			    >::value));
      while(beg < end)
	{
	  if (std::find_if(beg->second.begin(),beg->second.end(),
			   Sequence::invalidPolyChar())
	      != beg->second.end())
	    {
	      return false;
	    }
	  ++beg;
	}
      return true;
    }

    template < typename T >
    unsigned UnGappedLength (const std::vector <T>&data) 

    /*!
      Returns the number of sites in the alignment for which
      all objects do not contain the gap character '-'.  If the data are not 
      aligned, the value Sequence::SEQMAXUNSIGNED is returned as an error.
      \param data vector<T> to check
    */
    {
      BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			   || boost::is_same<std::pair<std::string,std::string>,T>::value));
      bool site_gapped = 0;
      int len = 0;
      if (!IsAlignment(data))
	return Sequence::SEQMAXUNSIGNED;

      for (size_t j = 0; j < data[0].second.length (); ++j)

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

    template <typename T>
    void RemoveGaps (std::vector <T> &data)
    /*!
      Modifies the data vector to remove all positions that
      contain the gap character'-'.
      \param data vector<T> to modify
    */
    {

       BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			    || boost::is_same<std::pair<std::string,std::string>,T>::value));
       size_t i, j,datasize=data.size();
       size_t length = data[0].second.length ();
      std::vector < std::string > ungapped_sequences(data.size());
      bool site_is_gapped;

      for (i = 0; i < length; ++i)
        {	//iterate over sites
          for ( j = 0, site_is_gapped = 0;
                j < datasize;  ++j)
            {
              if (data[j].second[i] == '-')
                {
                  site_is_gapped = 1;
                  j = datasize;
                }
            }
          if (!(site_is_gapped))
            {
              for ( j = 0 ; j != data.size();  ++j)
                ungapped_sequences[j] += data[j].second[i];
            }
        }

      //redo the data
      for (j = 0; j != datasize ; ++j)
        {
          data[j] = T(data[j].first,ungapped_sequences[j]);
        }
    }

    //only remove gaps from the beginning
    //and end of an alignment
    template < typename T >
    void RemoveTerminalGaps (std::vector <T>&data)
    /*!
      Remove all gapped sites from the ends of the alignment,
      up until the first site on either side that is ungapped.
      \param data vector<T> to modify
    */
    {
      BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			   || boost::is_same<std::pair<std::string,std::string>,T>::value));
      size_t i,j,length = data[0].second.length ();	//how much we have to iterate over
      std::vector < std::string > trimmed_sequences;
      size_t leftmost, rightmost, numUngapped;
      size_t offset;
      size_t size = data.size ();

      leftmost = SEQMAXUNSIGNED;
      rightmost = length + 1;	//offset by one b/c its an array...

      //find the leftmost site where all sites in the alignment are ungapped
      for (i = 0; i < length; ++i)
        {	//iterate over sites
          for (numUngapped = 0, j = 0; j != data.size (); ++j)
            {
              if (data[j].second[i] != '-')
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
      for (i = length - 1;  i < data[0].second.length() && exit_condition == false; --i)
        {
          for (numUngapped = 0, j = 0; j != data.size (); ++j)
            {
              if (data[j].second[i] != '-')
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
        trimmed_sequences.push_back (data[j].second.substr  (leftmost, rightmost-leftmost+1));

      //now, redo the seq array for the current object
      for ( j = 0; j != data.size (); ++j)
        {
          data[j] =T(data[j].first,trimmed_sequences[j]);
        }
    }

    template < typename T >
    void RemoveFixedOutgroupInsertions( std::vector<T> & data,
					unsigned site,
					const unsigned & ref )
    /*!
      Removes all positions from data that for which the outgroup
      contains an insertion relative to ingroup
      @param data a vector of Seq objects
      @param site index of the site at which to begin (set to 0 usually)
      @param ref the index of the outgroup in @a data
    */
    {
      BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			   || boost::is_same<std::pair<std::string,std::string>,T>::value));
      size_t nsam = data.size()-1;
      size_t s = data[0].second.length()-1,max=s+1;
      for(size_t i=0 ; i<max ; --s,++i)
	{
	  unsigned ngap=0;
	  for(size_t ind=0;ind<data.size();++ind)
	    {
	      if (ind != ref && data[ind].second[s] == '-')
		{
		  ngap++;
		}
	    }
	  if(ngap==nsam)
	    {
	      for(unsigned ind=0;ind<data.size();++ind)
		{
		  data[ind].second.erase(s,1);
		}
	    }
	}
    }

    template < typename T >
    std::vector < T >Trim (const std::vector < T >&data,
                           const std::vector <int> &sites) 

    /*!
      Returns a copy of the data vector, modified in the following way.  The sites vector
      contains an even number of sites (whose values are sorted). If sites 
      does not contain an even number of values Sequence::SeqException is thrown.
      If sites is empty, Sequence::SeqException is thrown.  The values in sites
      represent a series of intervals that you wish to keep, and the return vector is
      consists only of those--i.e. all positions not present
      in the intervals defined in sites are lost.  For example, if you pass a 
      vector<int> containing the values 0,10,21, and 30, then the data vector is modified
      so that positions 0 through 10 and 21 through 30 are all that remains.  
      One intended use of this function is to pull, for example, the coding region
      out of an aligned block.
      \param data the original data
      \param sites vector<int> containing an even number of integers 
      specifying the intervals of data to keep
      \exception Sequence::SeqException
    */
    {
      BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			   || boost::is_same<std::pair<std::string,std::string>,T>::value));
      size_t i, j, numseqs = data.size (), numIntervals =  sites.size ();
      unsigned start, stop;
      std::vector < T >trimmedData(numseqs);
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
              trimmedTemp[j] += data[j].second.substr (start, stop - start + 1);
            }
        }

      for (i = 0; i < numseqs; ++i)
        {
          trimmedData[i]=T(data[i].first,trimmedTemp[i]);
        }

      return trimmedData;
    }

    template < typename T >
    std::vector < T >TrimComplement (const std::vector < T >&data,
                                     const std::vector < int > &sites) 
      

    /*!
      Returns a copy the data vector, modified in the following way.  The sites vector
      contains an even number of sites (whose values are sorted). If sites 
      does not contain an even number of values Sequence::SeqException is thrown.
      If sites is empty, Sequence::SeqException is thrown.  The values in sites
      represent a series of intervals that you wish to keep, and the return vector 
      consists only of sites not present in \a sites--i.e. all positions not present
      in the intervals defined in sites are kept.  For example, if you pass a 
      vector<int> containing the values 0,10,21, and 30, then the data vector is modified
      so that positions 11 through 20 and 31 through the end of the sequences
      are all that remains.  
      \param data the original data
      \param sites vector<int> containing an even number of integers 
      specifying the intervals of data to throw away
      \exception Sequence::SeqException
    */
    {
      BOOST_STATIC_ASSERT((boost::is_base_and_derived<std::pair<std::string,std::string>,T>::value
			   || boost::is_same<std::pair<std::string,std::string>,T>::value));
      std::vector < int >newSites;
      size_t i, j, start, stop, numseqs = data.size (), numIntervals = sites.size (), lastval;

      if (sites.empty ())
        {
          throw SeqException ("Sequence::Alignment::TrimComplement(): empty vector of positions passed to function");
        }
      if (numIntervals % 2 != 0)
        {
          throw SeqException ("Sequence::Alignment::TrimComplement(): odd numer of positions passed to function");
        }

      std::vector < T >trimmedData(numseqs);
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
                data[j].second.substr (start,  stop - start + 1);
            }
        }

      for (j = 0; j < numseqs; ++j)
        {
          trimmedTemp[j] += data[j].second.substr (lastval);
        }

      for (i = 0; i < numseqs; ++i)
        {
          trimmedData[i] = T(data[i].first,trimmedTemp[i]);
        }

      return trimmedData;
    }

  }
}


