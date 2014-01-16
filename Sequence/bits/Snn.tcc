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

//  -*- C++ -*-
#ifndef __SEQUENCE_BITS_SNN_TCC__
#define __SEQUENCE_BITS_SNN_TCC__
#include <Sequence/Portability/random_shuffle.hpp>
#include <algorithm>

namespace Sequence
{
  template< typename uniform_int_generator >
  std::pair<double,double>
  Snn_test_details(const PolyTable & snpTable,
		   const unsigned config[],
		   const size_t & npop,
		   uniform_int_generator & uni_int,
		   const unsigned & nperms)
  {
    std::vector < std::vector<double> > dkj(snpTable.size(),
					    std::vector<double>(snpTable.size(),0.));
    for(unsigned i=0;i<snpTable.size()-1;++i)
      {
	for(unsigned j=i+1;j<snpTable.size();++j)
	  {
	    dkj[i][j] = Sequence::NumDiffs(snpTable[i],snpTable[j]);
	  }
      }
    std::vector <unsigned> __individuals(snpTable.size(),0u);
    for(unsigned i=0;i<snpTable.size();++i)
      {
	__individuals[i]=i;
      }
    const double observed = Snn_statistic(&__individuals[0],
					  dkj,
					  config,
					  npop,__individuals.size());
    unsigned perm = nperms;
    unsigned pv = 0;
    while( perm-- )
      {
#ifdef FORCE_STD_RANDOM_SHUFFLE
	std::random_shuffle(__individuals.begin(),__individuals.end(),uni_int);
#else
	Sequence::random_shuffle(__individuals.begin(),__individuals.end(),uni_int);
#endif
	double permuted_stat = Snn_statistic(&__individuals[0],
					     dkj,
					     config,
					     npop,__individuals.size());
	pv += (permuted_stat >= observed) ? 1 : 0;
      }
    return std::make_pair(observed,double(pv)/double(nperms));
  }

  template< typename uniform_int_generator >
  std::vector< std::vector<double> >
  Snn_test_pairwise_details(const PolyTable & snpTable,
			    const unsigned config[],
			    const size_t & npop,
			    uniform_int_generator uni_int,
			    const unsigned & nperms)
  {
    std::vector<unsigned> offsets(npop+1);
    for(unsigned i=1;i<=npop;++i)
      {
	offsets[i] = offsets[i-1]+config[i-1];
      }
    std::vector < std::vector<double> > dkj(snpTable.size(),
					    std::vector<double>(snpTable.size(),0.));
    for(unsigned i=0;i<snpTable.size()-1;++i)
      {
	for(unsigned j=i+1;j<snpTable.size();++j)
	  {
	    dkj[i][j] = Sequence::NumDiffs(snpTable[i],snpTable[j]);
	  }
      }
    std::vector <unsigned> __individuals(snpTable.size(),0u),individuals;
    unsigned _config[2];
    
    for(unsigned i=0;i<snpTable.size();++i)
      {
	__individuals[i]=i;
      }
    std::vector< std::vector<double> > rv;
    for(unsigned i=0;i< npop-1;++i)
      {
	for(unsigned j=i+1;j<npop;++j)
	  {
	    individuals.assign(__individuals.begin()+offsets[i],
			       __individuals.begin()+offsets[i+1]);
	    individuals.insert(individuals.end(),
			       __individuals.begin()+offsets[j],
			       __individuals.begin()+offsets[j+1]);
	    _config[0] = config[i];
	    _config[1] = config[j];
	    const double observed = Snn_statistic( &individuals[0],dkj,_config,2,
						   _config[0]+_config[1]); 
	    unsigned p =0, _nperms = nperms;
	    while(_nperms>0)
	      {
#ifdef FORCE_STD_RANDOM_SHUFFLE
		std::random_shuffle(__individuals.begin(),__individuals.end(),uni_int);
#else
		Sequence::random_shuffle(individuals.begin(),
					 individuals.end(),
					 uni_int);
#endif
		double perm = Snn_statistic( &individuals[0],dkj,_config,2,
					     _config[0]+_config[1]); 
		if(perm>=observed) ++p;
		_nperms--;
	      }
	    std::vector<double> result(4);
	    result[0]=double(i+1);
	    result[1]=double(j+1);
	    result[2]=double(observed);
	    result[3]=double(p)/double(nperms);
	    rv.push_back(result);
	  }
      }
    return rv;
  }

  template< typename uniform_int_generator >
  std::pair<double,double>
  Snn_test(const PolyTable & snpTable,
	   const unsigned config[],
	   const size_t & npop,
	   uniform_int_generator & uni_int,
	   const unsigned & nperms)
  /*!
    @brief Conducts a permutation-test of Hudson's Snn (sequence nearest-neighbor) statistic 
    \param snpTable The data on which we wish to perform the test.
    \param config An array of the sample sizes in each deme.
    \param npop The number of populations.  For example, npop could equal config.size() if config were a vector
    \param uni_int A random number generator whose operator() takes one argument, n, and returns a value uniformly-distributed on the half-open interval [0,N)
    \param nperms The number of permutations to do for the test
    \return A pair of doubles (std::pair<double,double>).  the first member of the pair is
    the observed value of the statistic, and the second member is the estimated p-value
    \ingroup popgenanalysis
  */
  {
    return Snn_test_details(snpTable,config,npop,uni_int,nperms);
  }

  template< typename uniform_int_generator >
  std::vector< std::vector<double> >
  Snn_test_pairwise(const PolyTable & snpTable,
		    const unsigned config[],
		    const size_t & npop,
		    uniform_int_generator & uni_int,
		    const unsigned & nperms)
  /*!
    @brief Conducts a permutation-test of Hudson's Snn (sequence nearest-neighbor) statistic,
    for all pairwise combinations of populations
    \param snpTable The data on which we wish to perform the test.
    \param config An array of the sample sizes in each deme.
    \param npop The number of populations.  For example, npop could equal config.size() if config were a vector
    \param uni_int A random number generator whose operator() takes one argument, n, 
    and returns a value uniformly-distributed on the half-open interval [0,N)
    \param nperms The number of permutations to do for the test
    \return A vector of vector<double>.  Each vector contains 4 elements, indexed 0 to 3.  Elements 0 and 1
    are the indexes of the i-th and j-th population.  Element 2 is the observed Snn between the i-th and j-th
    population, and element 3 is the p-value estimated by permutation
    \ingroup popgenanalysis
  */
  {
    return Snn_test_pairwise_details(snpTable,config,npop,uni_int,nperms);
  }

  template< typename uniform_int_generator >
  std::pair<double,double>
  Snn_test(const PolyTable & snpTable,
	   const unsigned config[],
	   const size_t & npop,
	   const uniform_int_generator & uni_int,
	   const unsigned & nperms)
  /*!
    @brief Conducts a permutation-test of Hudson's Snn (sequence nearest-neighbor) statistic 
    \param snpTable The data on which we wish to perform the test.
    \param config An array of the sample sizes in each deme.
    \param npop The number of populations.  For example, npop could equal config.size() if config were a vector
    \param uni_int A random number generator whose operator() takes one argument, n, and returns a value uniformly-distributed on the half-open interval [0,N)
    \param nperms The number of permutations to do for the test
    \return A pair of doubles (std::pair<double,double>).  the first member of the pair is
    the observed value of the statistic, and the second member is the estimated p-value
    \ingroup popgenanalysis
  */
  {
    return Snn_test_details(snpTable,config,npop,uni_int,nperms);
  }

  template< typename uniform_int_generator >
  std::vector< std::vector<double> >
  Snn_test_pairwise(const PolyTable & snpTable,
		    const unsigned config[],
		    const size_t & npop,
		    const uniform_int_generator & uni_int,
		    const unsigned & nperms)
  /*!
    @brief Conducts a permutation-test of Hudson's Snn (sequence nearest-neighbor) statistic,
    for all pairwise combinations of populations
    \param snpTable The data on which we wish to perform the test.
    \param config An array of the sample sizes in each deme.
    \param npop The number of populations.  For example, npop could equal config.size() if config were a vector
    \param uni_int A random number generator whose operator() takes one argument, n, 
    and returns a value uniformly-distributed on the half-open interval [0,N)
    \param nperms The number of permutations to do for the test
    \return A vector of vector<double>.  Each vector contains 4 elements, indexed 0 to 3.  Elements 0 and 1
    are the indexes of the i-th and j-th population.  Element 2 is the observed Snn between the i-th and j-th
    population, and element 3 is the p-value estimated by permutation
    \ingroup popgenanalysis
  */
  {
    return Snn_test_pairwise_details(snpTable,config,npop,uni_int,nperms);
  }
}//ns Sequence
#endif
