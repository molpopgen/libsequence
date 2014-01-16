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

#include <Sequence/FST.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/SeqConstants.hpp>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <limits>
/*!
  \class Sequence::FST Sequence/FST.hpp
  \ingroup popgenanalysis
  @short analysis of population structure using \f$F_{ST}\f$
*/

using std::accumulate;
using std::vector;
using std::pair;

namespace Sequence
{
#ifndef DOXYGEN_SKIP //doxygen should skip this

  struct FSTimpl
  {
    typedef std::vector<Sequence::stateCounter> _vCounts;
    typedef std::vector< _vCounts > _vvCounts;
    _vvCounts _Counts;
    mutable unsigned _nsam, _nsites, _npop;
    mutable double _piB, _piT, _piS, _piD;
    vector<unsigned> _config;
    vector<double> _weights;
    vector< pair<unsigned,unsigned> > _indexes;
    mutable bool _calcsDone;
    polySiteVector pv;
    typedef std::pair< std::set<char>, std::set<char> > PopStateSets;
    PopStateSets getPopStateSets(const unsigned &i,
				 const unsigned &j,
				 const unsigned & site);
    FSTimpl(const PolyTable *data, unsigned npop, const unsigned *config,
		     const double *weights, bool haveOutgroup, unsigned outgroup);
  };

  FSTimpl::FSTimpl(const PolyTable *data, unsigned npop, const unsigned *config,
	    const double *weights, bool haveOutgroup, unsigned outgroup) :
      _Counts( _vvCounts(npop,_vCounts(data->numsites()))),
      _nsam(unsigned(data->size())),
      _nsites(data->numsites()),_npop(npop),
      _piB(0.), _piT(0.), _piS(0.), _piD(0.),
      _indexes( vector< pair<unsigned,unsigned> >(npop) ),
      _calcsDone(false),
      pv( rotatePolyTable(data) )
    {
    if (config == NULL)
      {
        throw SeqException("Seqence::FST -- config vector is NULL");
      }
    else
      {
        _config.assign(config,config+npop);
        if(accumulate(_config.begin(),_config.end(),
		      unsigned(0))+unsigned(haveOutgroup) != _nsam)
          {
            throw SeqException("Seqence::FST -- sum of elements in config does not equal the sample size stored in data");
          }
      }
    if (weights == NULL)
      //equally weight all populations
      {
        _weights.assign(_npop, 1./double(_npop));
      }
    else
      {
        _weights.assign(weights,weights+npop);
	double sum = accumulate(_weights.begin(),_weights.end(),0.);
        if ( std::fabs(sum-1.) > std::numeric_limits<double>::epsilon() )
          {
            throw SeqException("Seqence::FST -- weights do not sum to 1");
          }
      }

    _indexes[0] = pair<unsigned,unsigned>(0,_config[0]);
    for(unsigned i = 1 ; i < _npop ; ++i)
      {
        _indexes[i]  =  pair<unsigned,unsigned>(_indexes[i-1].second,
						      _indexes[i-1].second+_config[i]);
      }
    for (unsigned i = 0 ; i < _npop ; ++i) //over pops
      {
        for(unsigned j = 0 ; j < _nsites ; ++j)//over sites
          {
            _Counts[i][j] = stateCounter('-');
            for(unsigned k = _indexes[i].first ; k < _indexes[i].second ; ++k)//over sequences
              {
                if (haveOutgroup == false ||
                    (haveOutgroup == true && k != outgroup))
		  _Counts[i][j]((*data)[k][j]);
              }
          }
      }

    if (haveOutgroup == true)
      {
	//remove outgroup state from pv
	for(unsigned i=0;i<pv.size();++i)
	  {
	    pv[i].second.erase(outgroup,1);
	  }
      }
    }

  FSTimpl::PopStateSets FSTimpl::getPopStateSets(const unsigned &pop1,
						 const unsigned &pop2,
						 const unsigned & site)
  {
    size_t beg1 = std::accumulate(_config.begin(),
				  _config.begin()+pop1,
				  0u);
    size_t size1 = *(_config.begin()+pop1);
    size_t beg2 = std::accumulate(_config.begin(),
				  _config.begin()+pop2,
				  0u);
    size_t size2 = *(_config.begin()+pop2);
    
    //we need some rigamarole here to skip missing data
    vector<char> ch1(pv[site].second.begin()+beg1,
		     pv[site].second.begin()+beg1+size1),
      ch2(pv[site].second.begin()+beg2,
	  pv[site].second.begin()+beg2+size2);
    ch1.erase( std::remove(ch1.begin(),ch1.end(),'N'),ch1.end() );
    ch2.erase( std::remove(ch2.begin(),ch2.end(),'N'),ch2.end() );
    std::set<char> states1(ch1.begin(),ch1.end()),
      states2(ch2.begin(),ch2.end());
    return std::make_pair(states1,states2);
  }
#endif
  FST::FST(const PolyTable *data, unsigned npop, const unsigned *config,
           const double *weights, bool haveOutgroup, unsigned outgroup) 


    /*!
      \param data a Sequence::PolyTable with the data for all populations
      \param npop the number of populations
      \param config a list of sample sizes.  the number of elements in \a config must equal \a npop, 
      otherwise a segfault may occur.  An exception is thrown if the sum of the elements in config
      does not equal \a data->size(). \note If \a haveOutgroup is true, then the sum of the elements in 
      \a config must be 1 less than the size of \a data (i.e. it should sum to the size of the 
      polymorphism sample only!).
      \param weights a list of weights to assign to each population.  If NULL, all pops are weighted
      equally.  The number of elements in weights must be >= \a npop.  Also, the sum of the elements
      in the range &weights[0] to &weights[npop] must be 1, else an exception will be thrown.
      \param haveOutgroup false if no outgroup sequence is present, set to true otherwise. Please note
      that it is rather unwise to have an outgroup tested when doing a permutation test of Fst, since
      the index of the outgroup will likely get scrambled and not correspond to \a outgroup.
      \param outgroup if \a haveOutgroup is true, \a outgroup is the index of the outgroup in \a data
    */
  {
    try
      {
	impl = std::auto_ptr<FSTimpl>(new FSTimpl(data,npop,config,weights,haveOutgroup,outgroup));
      }
    catch(Sequence::SeqException &e)
      { 
	throw;
      }
    catch (...)
      {
	throw (Sequence::SeqException("Sequence::FST : unknown exception caught by constructor"));
      }
  }

  FST::~FST(void)
  {
    //delete impl;
  }

  void FST::doCalcs(void) const
  {
    //calculate sums of w/in population weights and heterozygosity
    double w_ii_sq = 0., weighted_Pi_ii = 0.;

    for (unsigned i = 0 ; i < impl->_npop ; ++i) //over pops
      {
	double Pi = 0.;
	w_ii_sq += impl->_weights[i]*impl->_weights[i];
	for(unsigned j = 0 ; j < impl->_nsites ; ++j)//over sites
	  {
	    unsigned n = impl->_config[i];
	    double SSH = 0.;
	    n -= impl->_Counts[i][j].n;
	    double denom = (double(n)* (double(n) - 1.0));
	    SSH += (impl->_Counts[i][j].a > 0) ? double(impl->_Counts[i][j].a) * 
	      double (impl->_Counts[i][j].a-1) /denom : 0. ;
	    SSH += (impl->_Counts[i][j].g > 0) ? double(impl->_Counts[i][j].g) * 
	      double (impl->_Counts[i][j].g-1) /denom : 0. ;
	    SSH += (impl->_Counts[i][j].c > 0) ? double(impl->_Counts[i][j].c) * 
	      double (impl->_Counts[i][j].c-1) /denom : 0. ;
	    SSH += (impl->_Counts[i][j].t > 0) ? double(impl->_Counts[i][j].t) * 
	      double (impl->_Counts[i][j].t-1) /denom : 0. ;
	    SSH += (impl->_Counts[i][j].zero > 0) ? double(impl->_Counts[i][j].zero) *
	      double (impl->_Counts[i][j].zero-1) /denom : 0. ;
	    SSH += (impl->_Counts[i][j].one > 0) ? double(impl->_Counts[i][j].one) *
	      double (impl->_Counts[i][j].one-1) /denom : 0. ;
	    Pi += (1.0 - SSH);
	  }//sites
	weighted_Pi_ii += impl->_weights[i]*impl->_weights[i]*Pi;
      }//pops

    //now calculate between-population divergence,
    //using stored state info in impl->_Counts
    double sum_wi_wj = 0., weighted_Pi_ij = 0.;
    for(unsigned i = 0 ; i < impl->_npop - 1 ; ++i)//over pops_i
      {
	for(unsigned j = i+1 ; j < impl->_npop ; ++j)//over pops_j
	  {
	    double Pi_i_j = 0.;
	    sum_wi_wj += impl->_weights[i]*impl->_weights[j];
	    for(unsigned k = 0 ; k < impl->_nsites ; ++k)//over sites
	      {
		unsigned ni = impl->_config[i],nj = impl->_config[j];
		ni -= impl->_Counts[i][k].n;
		nj -= impl->_Counts[j][k].n;
		Pi_i_j += (impl->_Counts[i][k].a > 0) ?
		  (double(impl->_Counts[i][k].a)/double(ni))*
		  (double(nj-impl->_Counts[j][k].a)/double(nj)):0.;
		Pi_i_j += (impl->_Counts[i][k].g > 0) ?
		  (double(impl->_Counts[i][k].g)/double(ni))*
		  (double(nj-impl->_Counts[j][k].g)/double(nj)):0.;
		Pi_i_j += (impl->_Counts[i][k].c > 0) ?
		  (double(impl->_Counts[i][k].c)/double(ni))*
		  (double(nj-impl->_Counts[j][k].c)/double(nj)):0.;
		Pi_i_j += (impl->_Counts[i][k].t > 0) ?
		  (double(impl->_Counts[i][k].t)/double(ni))*
		  (double(nj-impl->_Counts[j][k].t)/double(nj)):0.;
		Pi_i_j += (impl->_Counts[i][k].zero > 0) ?
		  (double(impl->_Counts[i][k].zero)/double(ni))*
		  (double(nj-impl->_Counts[j][k].zero)/double(nj)):0.;
		Pi_i_j += (impl->_Counts[i][k].one > 0) ?
		  (double(impl->_Counts[i][k].one)/double(ni))*
		  (double(nj-impl->_Counts[j][k].one)/double(nj)):0.;
	      }//sites
	    weighted_Pi_ij += impl->_weights[i]*impl->_weights[j]*Pi_i_j;
	  }//pops_j
      }//pops_i

    //calculate measures of diversity/divergence needed to obtain Fst
    impl->_piT = weighted_Pi_ii + 2.*weighted_Pi_ij;
    impl->_piS = weighted_Pi_ii / w_ii_sq ;
    impl->_piB = weighted_Pi_ij / sum_wi_wj;
    impl->_piD = (impl->_piT - impl->_piS)/(2. * sum_wi_wj);

    impl->_calcsDone = true;
  }

  std::set<double> FST::shared(unsigned pop1, unsigned pop2) const
    /*!
      \return an object of type std::set<double> representing the positions
      where shared polymorphisms are inferred
      \param pop1 a population index
      \param pop2 a population index
    */
  {
    if (pop1 > impl->_npop-1 || pop2 > impl->_npop-1)
      return std::set<double>();

    std::set<double> sharedList;
    
    for(unsigned site = 0 ; site < impl->_nsites ; ++site)
      {
	if (impl->_Counts[pop1][site].gap == 0 &&
	    impl->_Counts[pop2][site].gap == 0)
	  {
	    FSTimpl::PopStateSets states = impl->getPopStateSets(pop1,pop2,site);
	    if (states.first.size() > 1 && states.second.size() > 1)
	      {
		vector<char> overlap(states.first.size()+states.second.size());
		vector<char>::iterator itr = std::set_intersection(states.first.begin(),
								   states.first.end(),
								   states.second.begin(),
								   states.second.end(),
								   overlap.begin());
		if ( itr - overlap.begin() > 0) //shared states exist
		  {
		    sharedList.insert(impl->pv[site].first);
		  }
	      }
	  }
      }
    return sharedList;
  }

  std::set<double> FST::fixed(unsigned pop1, unsigned pop2)  const
    /*!
      \return an object of type std::set<double> representing the positions of sites at which
      fixed differences occur and a list of positions at which fixed differences occur
      \param pop1 a population index
      \param pop2 a population index
    */
  {
    if (pop1 > impl->_npop-1 || pop2 > impl->_npop-1)
      return std::set<double>();
    std::set<double> fixedList;
    for(unsigned site = 0 ; site < impl->_nsites ; ++site)
      {
	if (impl->_Counts[pop1][site].gap == 0 &&
	    impl->_Counts[pop2][site].gap == 0)
	  {
	    FSTimpl::PopStateSets states = impl->getPopStateSets(pop1,pop2,site);
	    vector<char> overlap(states.first.size()+states.second.size());
	    vector<char>::iterator itr = std::set_intersection(states.first.begin(),
							       states.first.end(),
							       states.second.begin(),
							       states.second.end(),
							       overlap.begin());
	    if ( itr - overlap.begin() == 0) //no shared states exist
	      {
		fixedList.insert(impl->pv[site].first);
	      }
	  }
      }
    return fixedList;
  }

  std::pair< std::set<double>,std::set<double> >
  FST::Private(unsigned pop1, unsigned pop2) const
    /*!
      \return std::pair containing two objects of type std::set<double>
      The first member of the pair represents the segregating sites with 
      private mutations in \a pop1,
      the second member the sites with private polymorphisms in \a pop2. 
      Returns a pair of empty sets if \a pop1  or \a pop2 are out of range
      \param pop1 a population index
      \param pop2 a population index
    */
  {
    if (pop1 > impl->_npop-1 || pop2 > impl->_npop-1)
      return std::pair< std::set<double>,std::set<double> >();

    std::set<double> p1,p2;
    for(unsigned site = 0 ; site < impl->_nsites ; ++site)
      {
	if (impl->_Counts[pop1][site].gap == 0 &&
	    impl->_Counts[pop2][site].gap == 0)
	  {
	    FSTimpl::PopStateSets states = impl->getPopStateSets(pop1,pop2,site);

	    //now, use std::set_difference to identify private polys
	    vector<char> priv1(states.first.size()),priv2(states.second.size());
	    //private in pop1
	    vector<char>::iterator itr1 = std::set_difference(states.first.begin(),states.first.end(),
							      states.second.begin(),states.second.end(),
							      priv1.begin());
	    //private in pop2
	    vector<char>::iterator itr2 = std::set_difference(states.second.begin(),states.second.end(),
							      states.first.begin(),states.first.end(),
							      priv2.begin());
	    //if there are unique states && site is polymorphic,
	    //then there is a private poly
	    if ( (itr1-priv1.begin())>0 && states.first.size()>1)
	      {
		p1.insert(impl->pv[site].first);
	      }
	    if ( (itr2-priv2.begin())>0 && states.second.size()>1)
	      {
		p2.insert(impl->pv[site].first);
	      }
	  }
      }
    return std::make_pair(p1,p2);
  }

  double FST::HSM(void) const
    /*!
      \return \f[F_{ST}=\frac{\pi_{D}}{\pi_S + \pi_D},\f] which is the
      definition of \f$F_{ST}\f$ according to Hudson, Slatkin and Maddison (1992)
      Estimation of levels of gene flow from population data. Genetics 132:583-589
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return impl->_piD/(impl->_piS+impl->_piD);
  }

  double FST::Slatkin(void) const
    /*!
      \return \f[F_{ST}=\frac{\pi_D}{2\pi_S + \pi_D},\f] which is the
      definition of \f$F_{ST}\f$ according to Slatkin (1993) Isolation by distance in
      equilibrium and non-equilibrium populations. Evolution 47: 264-279
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return impl->_piD/(2.*impl->_piS + impl->_piD);
  }

  double FST::HBK(void) const
    /*!
      \return \f[F_{ST}= 1 - \frac{\pi_S}{\pi_T} , \f] which is the
      definition of \f$F_{ST}\f$ according to Hudson, Boos, and Kaplan (1992)
      A statistical test for detecting geographic subdivision. Mol. Biol. Evol.
      9:138-151
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return 1.-(impl->_piS/impl->_piT);
  }

  double FST::piB(void) const
    /*!
      \return \f[\pi_B= \frac{\sum_{i<j}w_i w_j \pi_{ij}}{\sum_{i<j}w_i w_j}, \f]
      which is the mean parwise divergence between 2 alleles drawn from 2 populations
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return impl->_piB;
  }

  double FST::piT(void) const
    /*!
      \return \f[\pi_T = \sum_i w_i^2 \pi_{ii} + 2\sum_{i<j}w_i w_j \pi_{ij}, \f] which is the 
      total diversity in the sample
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return impl->_piT;
  }

  double FST::piS(void) const
    /*!
      \return \f[\pi_S = \frac{\sum_i w_i^2 \pi_{ii}}{\sum_i w_i^2},\f] which is the mean within-population
      diversity
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return impl->_piS;
  }

  double FST::piD(void) const
    /*!
      \return \f[ \pi_D = \frac{\pi_T - \pi_S}{2 \sum_{i<j}w_i w_j}, \f] which is a measure of
      between-population divergence that is proportional to \f$t_1 - t_0,\f$ where \f$t_1\f$ is the
      mean time to coalescence for 2 alleles drawn from different populations, and \f$t_o\f$ is the mean time
      to coalescence for 2 alleles drawn from the same population
    */
  {
    if(impl->_calcsDone == false)
      doCalcs();
    return impl->_piD;
  }
}
