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

/*! \file PolyFunctional.hpp
  @brief A collection of function objects for SNP analysis.

  These objects may be useful for certain one-off calculations
  in lieu of Sequence::PolySNP.  They are also useful guides
  for how to implement functors that can be used in conjunction
  with STL algorithms to do calculations on objects of type
  Sequence::PolyTable.  This file declares Sequence::countStates,
  Sequence::countDerivedStates,Sequence::ssh, and Sequence::nmuts.
*/
#ifndef __POLY_FUNCTIONAL_HPP__
#define __POLY_FUNCTIONAL_HPP__

#include <Sequence/stateCounter.hpp>
#include <Sequence/PolyTableManip.hpp>
#include <Sequence/typedefs.hpp>
#include <vector>
#include <string>
#include <utility>

namespace Sequence
{
  struct countStates
  /*!
    \brief Functor to count the number of states, excluding gaps and missing data,
    in a range of characters.
  */
  {
    template<typename charItr>
    inline stateCounter operator()( charItr beg,
				    charItr end,
				    const bool & haveOutgroup = false,
				    const size_t & outgroup = 0 )
    /*!
      \param beg beginning of a range of characters
      \param end end of a range of characters
      \param haveOutgroup true of one of the elements of \a site is an outgroup state,
      false otherwise
      \param outgroup the index of the outgroup sequence in \a site
      \note will return an empty object if beg>=end or beg+outgroup>=end
    */
    {
      stateCounter c;
      if(beg>=end || beg+outgroup>=end) return c;
      for( charItr itr=beg ; 
	   itr != end ; 
	   ++itr )
	{
	  if ( (!haveOutgroup) 
	       //casting to unsigned is safe since itr is always >= beg
	       || ( haveOutgroup && ( unsigned(itr-beg) != outgroup ) ) )
	    {
	      c(*itr);
	    }
	}
      return c;
    }
  };
  
  struct countDerivedStates
  /*!
    \brief Functor to count the number of derived states, excluding gaps and missing data,
    in a range of characters.
  */
  {
    template<typename charItr>
    inline stateCounter operator()( charItr beg,
				    charItr end,
				    const bool & haveOutgroup = true,
				    const size_t & outgroup = 0 )
    /*!
      \param beg beginning of a range of characters
      \param end end of a range of characters
      \param haveOutgroup true of one of the elements of \a site is an outgroup state,
      false otherwise
      \param outgroup the index of the outgroup sequence in \a site
      \note will return an empty object if beg>=end or beg+outgroup>=end or haveOutgroup==false
    */
    {
      stateCounter c;
      if(beg>=end || beg+outgroup>=end || haveOutgroup == false) return c;
      for( charItr itr=beg ; 
	   itr != end ; 
	   ++itr )
	{
	  //casting to unsigned is safe since itr is always >= beg
	  if ( (unsigned(itr-beg) != outgroup) 
	       && (*itr != *(beg+outgroup))  )
	    {
	      c(*itr);
	    }
	}
      return c;
    }
  };

  /*! \struct ssh Sequence/PolyFunctional.hpp
    A function object to keep track of "ssh", the sum of 
    site heterozygosity, aka pi, aka nucleotide diversity.
    It can be used with std::accumulate to do calculations
    accross a range of sites, i.e.:
    \code
    #include <numeric>
    #include <algorithm>
    #include <iostream>
    #include <boost/bind.hpp>
    #include <Sequence/FastaExplicit.hpp>
    #include <Sequence/PolySites.hpp>
    #include <Sequence/PolySNP.hpp>
    #include <Sequence/stateCounter.hpp>
    #include <Sequence/PolyFunctional.hpp>
		 
    int main()
    {
    std::vector< Sequence::Fasta > data;
    Sequence::Alignment::GetData(data,std::cin);
    Sequence::PolySites p(data);
    //calculate nucleotide diversity using STL
    double pi = std::accumulate(p.sbegin(),p.send(),0.,ssh());
    std::cout << pi
    << '\n';
    //calculate nucleotide diversity using Sequence::PolySNP
    Sequence::PolySNP a(&p);
    std::cout << a.ThetaPi() 
    << '\n';
    }
    \endcode
    \brief Calculate nucleotide diversity from a polymorphic site
  */
  struct ssh : public std::binary_function< double,polymorphicSite,double >

  {
    inline double operator()( double & _ssh,
			      const polymorphicSite & site,
			      const bool & haveOutgroup = false,
			      const unsigned & outgroup = 0 ) const
      /*!
	\param _ssh a value of ssh to increment
	\param site an object representing the value type of 
	PolyTable::const_site_iterator
	\param haveOutgroup true of one of the elements of \a site is an outgroup state,
	false otherwise
	\param outgroup the index of the outgroup sequence in \a site
      */
      {
	stateCounter c = countStates()(site.second.begin(),
				       site.second.end(),
				       haveOutgroup,outgroup);

	unsigned nsam = site.second.length() - 
	  unsigned( haveOutgroup==true ? 1 : 0 ) - c.n;
	double hom=0.;
	if (c.gap==0 && nsam > 1)
	  {
	    hom += (c.a > 0) ? double(c.a) * double(c.a-1) : 0.;
	    hom += (c.g > 0) ? double(c.g) * double(c.g-1) : 0.;
	    hom += (c.c > 0) ? double(c.c) * double(c.c-1) : 0.;
	    hom += (c.t > 0) ? double(c.t) * double(c.t-1) : 0.;
	    hom += (c.zero > 0) ? double(c.zero) * double(c.zero-1) : 0.;
	    hom += (c.one > 0) ? double(c.one) * double(c.one-1) : 0.;
	  }
	return _ssh += (hom > 0.) ? 1.- (hom/( double(nsam)*double(nsam-1) )) : 0.;
      }
  };

  /*! \class nmuts Sequence/PolyFunctional.hpp
    Function object to keep track of the number of mutations
    in a range of polymorphic sites. The template argument
    specifies the counting behavior, i.e. countStates
    or countDerivedStates.  It is designed to be used with
    std::accumulate, and will need boost::bind for the necessary arguments
    \code
    #include <numeric>
    #include <algorithm>
    #include <iostream>
    #include <boost/bind.hpp>
    #include <Sequence/FastaExplicit.hpp>
    #include <Sequence/PolySites.hpp>
    #include <Sequence/PolySNP.hpp>
    #include <Sequence/stateCounter.hpp>
    #include <Sequence/PolyFunctional.hpp>
		 
    int main()
    {
    std::vector< Sequence::Fasta > data;
    Sequence::Alignment::GetData(data,std::cin);
    Sequence::PolySites p(data);
    //calculate # mutations via algorithm
    unsigned nm = std::accumulate(p.sbegin(),p.send(),0u,
    boost::bind(nmuts<countStates>(),_1,_2,
    false,0));
    std::cout << nm
    << '\n';
    //for comparison, use Sequence::PolySNP
    Sequence::PolySNP a(&p);
    std::cout << a.NumMutations() << '\n';
    }
    \endcode
    \brief Calculate the number of mutations at a polymorphic site
  */
  template<typename counter>
  struct nmuts 
  {
    /*!
      allows boost::bind to be used in a simple way
    */
    typedef unsigned result_type;
    inline unsigned operator()(unsigned & nm,
			       const polymorphicSite & site,
			       const bool & haveOutgroup = false,
			       const unsigned & outgroup = 0 ) const
    /*!
      \param nm a value of nmuts to increment
      \param site an object representing the value type of 
      PolyTable::const_site_iterator
      \param haveOutgroup true of one of the elements of \a site is an outgroup state,
      false otherwise
      \param outgroup the index of the outgroup sequence in \a site
    */
    {
      stateCounter c = counter()(site.second.begin(),
				 site.second.end(),
				 haveOutgroup,outgroup);
      unsigned n = 0;
      n += (c.a > 0) ? 1 : 0;
      n += (c.g > 0) ? 1 : 0;
      n += (c.c > 0) ? 1 : 0;
      n += (c.t > 0) ? 1 : 0;
      n += (c.zero > 0) ? 1 : 0;
      n += (c.one > 0) ? 1 : 0;
      return nm += (n>=2) ? n-1 : 0;
    }
  };
}
#endif
