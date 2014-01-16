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

#ifndef __HKA_HPP__
#define __HKA_HPP__

#include <vector>
#include <boost/tuple/tuple.hpp>

namespace Sequence
{
  /*! \struct HKAdata Sequence/HKA.hpp
    Data from a single locus for an HKA test
    \short Data from a single locus for an HKA test
    \ingroup popgen
  */
  struct HKAdata
  {
    /*!
      SA is the number of polymorphic sites in species A,
      SB is the same quantity for species B
    */
    mutable unsigned SA,SB; //# seg sites
    /*!
      D is a measure of divergence between species A and B,
      i.e. mean pairwise differences
    */
    mutable double D;       // divergence
    /*!
      The sample sizes in species A and B
    */
    mutable unsigned nA,nB; //sample sizes
    HKAdata();
    HKAdata(const HKAdata & d);
    HKAdata(unsigned sa,unsigned sb,double d,
	    unsigned na,unsigned nb);
  };

  /*! \struct HKAresults Sequence/HKA.hpp
    Store the list of estimates of theta for each locus,
    as well as f (ratio of Ne of spp B relative to sppA), and T,
    the divergence time between species (on the scale of the
    Ne of spp A).
    \short results of calculations of the HKA test
    \ingroup popgen
  */
  struct HKAresults
  {
    /*!
      Vector of theta estimates from the HKA test. They are in the same
      order as the loci that were placed in the std::vector< HKAdata >
      passed to calcHKA
    */
    const std::vector< double > thetas;
    /*!
      The tuple represents the contribution of a locus
      to the total value of the HKA chi-squared statistic
    */
    typedef boost::tuple<double,double,double,double> chisq_tuple;
    /*!
      The enum documents the order in which the deviations
      are stored in chisq_tuple
    */
    enum chisq_tuple_elements {POLY,DIV,POLYA,POLYB};
    /*!
      The order of elements is the same as the loci that were 
      placed in the std::vector< HKAdata > passed to calcHKA.
      Members of the tuple are accessed via the boost::get<int>
      template function.  For example, if you have an HKAresults
      object called hkares, and you want to retrieve the deviation
      in divergence of the first locus:
      \code
      boost::get< DIV >(hkares.chisquareds[0]);
      \endcode
      The template argument to boost::get can be taken from 
      chisq_tuple_elements so that it's readable.
    */
    std::vector< chisq_tuple > chisquareds;
    /*!
      The relative Ne of species B relative to
      species A
    */
    const double fhat;
    /*!
      The estimate of the species divergence time,
      on the scale (1+fhat)*4Ne generations
    */
    const double That;
    /*!
      The HKA chi-squared statistic
    */
    const double xsq;
    /*!
      The HKA chi-squared statistic,
      only considering species A
    */
    const double xsqA;
    /*!
      The HKA chi-squared statistic,
      only considering species B
    */
    const double xsqB;
    explicit HKAresults (const std::vector<double> & _thetas,
			 const std::vector< chisq_tuple > & _chisquareds,
			 const double & _fhat, const double & _That,
			 const double & _xsq,
			 const double & _xsqA, const double & _xsqB);
  };

  HKAresults calcHKA ( const std::vector< HKAdata > & data );
} //namespace Sequence

#endif
