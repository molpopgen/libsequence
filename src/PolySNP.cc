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

#include <cmath>
#include <cfloat>
#include <cassert> 
#include <utility>
#include <cstdlib>
#include <cctype>
#include <set>
#include <algorithm>
#include <boost/bind.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/Comparisons.hpp>
#include <Sequence/Recombination.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolySNPimpl.hpp>
#include <Sequence/StringComp.hpp>

/*!
  \defgroup popgenanalysis Analysis of molecular population genetic data
  \ingroup popgen
*/

using std::string;
using namespace Sequence::Recombination;

namespace Sequence
{
  struct uniqueSeq : public std::binary_function<std::string,std::string,bool>
 {
    inline bool operator()(const std::string & l, const std::string  & r) const
    {
      //use Sequence::Different to prevent missing sites 
      //causing 2 sequences to be labelled as distinct
      return (  Different(l,r) && std::lexicographical_compare(l.begin(),l.end(),r.begin(),r.end(),lt_nocase()) );
    }
  };
  _PolySNPImpl::_PolySNPImpl (const Sequence::PolyTable * data, bool haveOutgroup ,
			      unsigned outgroup, bool totMuts):
    _data(data),
    _nsites(data->numsites()),
    _nsam(unsigned(data->size())),
    _outgroup(outgroup),
    _haveOutgroup(haveOutgroup),
    _totMuts(totMuts),
    _totsam (unsigned(data->size())),
    _DVK(0),
    _DVH(1.0),
    _counted_singletons(false),
    _know_pi(false),
    _CalculatedDandV(false),
    _pi(0.),
    _singletons(0),
    _walls_Bprime(0),
    _NumPoly(0),
    _walls_B(0.),_walls_Q(0.),
    _calculated_wall_stats(false),
    _counts(std::vector< Sequence::stateCounter >(_nsites,Sequence::stateCounter('-'))),
    _derivedCounts(std::vector< std::pair< bool,  stateCounter > >
		   (_nsites,
		    std::make_pair<bool,stateCounter>(true,stateCounter('-')))),
    _preprocessed(false)
  {
    if (haveOutgroup)
      --_totsam;	//because one sequence in data is an outgroup!
    preprocess();
  }

  void  _PolySNPImpl::preprocess()
  /*!
    This routine takes the data and obtains count information
    for each possible character state at each site.
    The reason for doing this is that the various summary statistics
    that depend on frequency spectrum information are all O(n x S) to
    compute, where n is the sample size and S the number of sites.
    However, all of these statistics depend basically on the same
    features of the data, so it inefficient to recalculate them
    every time.  This function obtains the counts for each site.
    Note that the efficiency gain is huge (especially for large
    data sets) as each summary statistic is reduced from an O(n x S)
    to an O(S) calculation.
    The increase in run-time efficiency comes at the cost
    of allocating 2 vectors whose sizes are linear in S
  */
  {
    if (!_preprocessed)
      {
	for(unsigned site = 0 ; site < _nsites ; ++site)
	  {
	    for (unsigned seq = 0 ; seq < _nsam ; ++seq)
	      {
		//process counts w/o respect to
		//ancestral or derived
		if (!_haveOutgroup
		    || (_haveOutgroup && seq != _outgroup))
		  {
		    _counts[site]( (*_data)[seq][site] );
		  }
		//process derived states if outgroup is
		//present
		if (_haveOutgroup == true)
		  {
		    //if outgroup state is missing data
		    //or a gap
		    //set the bool for that site to false
		    if ( std::toupper((*_data)[_outgroup][site]) == 'N' ||
			 (*_data)[_outgroup][site] == '-')
		      {
			_derivedCounts[site].first = false;
		      }
		    else
		      {
			//tabulate the derived state
			_derivedCounts[site].first = true;
			if ( seq != _outgroup && (*_data)[seq][site] != (*_data)[_outgroup][site] )
			  {
			    _derivedCounts[site].second((*_data)[seq][site]);
			  }
		      }
		  }
		else
		  {
		    //if no outgroup, set bool
		    //for site to false
		    _derivedCounts[site].first = false;
		  }
	      }
	    if(_counts[site].nStates() > 1 && _counts[site].gap==0) ++_NumPoly;
	  }
	_preprocessed = true;
      }
  }

  PolySNP::PolySNP (const Sequence::PolyTable * data, bool haveOutgroup ,
                    unsigned outgroup, bool totMuts):
    rep(std::auto_ptr<_PolySNPImpl>(new _PolySNPImpl(data,haveOutgroup,outgroup,totMuts)))
    /*!
      \param data a valid object of type Sequence::PolyTable
      \param haveOutgroup \c true if an outgroup is present, \c false otherwise
      \param outgroup if \a haveOutgroup is \c true, \a outgroup is the index of that sequence in data
      \param totMuts if true (the default) use the total number of inferred mutations, 
      otherwise use the total number of polymorphic sites in calculations
      \note this constructor allocates a pointer to an implementation class, which 
      automatically pre-processes the SNP data
    */
  {}

  PolySNP::~PolySNP (void)
  {
  }

  double
  PolySNP::ThetaPi (void)
  /*!
    Calculated here as the sum of 1.0 - sum of site homozygosity accross sites.\n 
    \f[
    \widehat\theta_\pi=\sum_{i=1}^{i=S}(1-\sum_{j=1}^{j=4}\frac{k_{j,i} \times (k_{j,i}-1)}{n_i \times (n_i - 1)});k_{j,i}>0
    \f]
    Where \f$S\f$ is the number of segregating sites, \f$k_{j,i}\f$ is the number of occurences of the \f$j_{th}\f$
    character state at site \f$i\f$, and \f$n_i\f$ is the sample size at site \f$i\f$.
    Calculating the statistic in this manner makes it easy to generalize to an arbitrary number
    of character states per polymorphic site\n
    Also equivalent to sum of site heterozygosities:\n
    \f[
    \widehat\theta_\pi=\sum_{i=1}^{i=S}2{p_i}{q_i}
    \f]\n
    Also equivalent to mean pairwise differences, but that's slow to calculate.\n
    \n
    If there is missing data (indicated by 'N' characters), the sample size
    is reduced for that site.  For example, if the data for the \f$i_{th}\f$ site is:\n
    A\n
    A\n
    A\n
    N\n
    N\n
    G\n
    Then ThetaPi is calculated for that site as if the sample size were 4 (not 6),
    and the polymorphic site frequencies are 3/4 for A and 1/4 for G\n\n
  */
  {
    assert ( rep->_preprocessed );
    if (rep->_know_pi == false)
      {
        double Pi = 0.0;
        for (unsigned i = 0;  i < rep->_nsites; ++i)
          {			//iterate over sites
	    if ( rep->_counts[i].gap == 0 &&
		 rep->_counts[i].nStates() > 1 )
	      {
		unsigned samplesize = rep->_totsam;
		double SSH = 0.0;	//sum of site homozygosity
		samplesize -= rep->_counts[i].n; //adjust sample size for missing data
		if (samplesize > 1)
		  {
		    double denom = (double(samplesize)* (double(samplesize) - 1.0));
		    SSH += (rep->_counts[i].a > 0) ? double(rep->_counts[i].a) * 
		      double (rep->_counts[i].a-1) /denom : 0. ;
		    SSH += (rep->_counts[i].g > 0) ? double(rep->_counts[i].g) * 
		      double (rep->_counts[i].g-1) /denom : 0. ;
		    SSH += (rep->_counts[i].c > 0) ? double(rep->_counts[i].c) * 
		      double (rep->_counts[i].c-1) /denom : 0. ;
		    SSH += (rep->_counts[i].t > 0) ? double(rep->_counts[i].t) * 
		      double (rep->_counts[i].t-1) /denom : 0. ;
		    SSH += (rep->_counts[i].zero > 0) ? double(rep->_counts[i].zero) * 
		      double (rep->_counts[i].zero-1) /denom : 0. ;
		    SSH += (rep->_counts[i].one > 0) ? double(rep->_counts[i].one) * 
		      double (rep->_counts[i].one-1) /denom : 0. ;
		    Pi += (1.0 - SSH);
		  }
	      }
          }
        rep->_pi = Pi;
        rep->_know_pi = true;
        return rep->_pi;
      }
    else
      return rep->_pi;
    return rep->_pi;
  }


  double
  PolySNP::ThetaW (void)
  /*!
    The classic "Watterson's Theta" statistic, generalized to missing data
    and multiple mutations per site:
    \f[
    \widehat\theta_w=\sum_{i=1}^{i=S}\frac{S}{\sum_{j=1}^{j=n_i-1}\frac{1}{j}}
    \f]\n
    For this statistic, \f$S\f$ is either the number of segregating sites,
    or the number of mutations on the genealogy and \f$n_i\f$ is the sample size
    at site i. If totMuts == 1,
    the number of mutations is used, else the number ofsegregating sites is used.
    \warning Statistic undefined if there are untyped SNPs.  In the presence of
    missing data, ThetaW is calculated as the sum (over all segregating sites)
    of 1/a_sub_n, where a_sub_n is the denominator of ThetaW, using the number 
    of alleles without missing data as the sample size.  More formally, the
    routine returns a calculation base on an unweighted sample size adjustment accross
    sites.
  */
  {
    assert ( rep->_preprocessed );
    double W = 0.0;
    for (unsigned i = 0;  i < rep->_nsites; ++i)
      {			//iterate over sitesvv
	if ( rep->_counts[i].gap == 0)
	  {
	    int nstates = rep->_counts[i].nStates();
	    unsigned nsam_site = rep->_totsam - rep->_counts[i].n;
	    double denom=0.0;
	    if(rep->_totMuts==true&&nstates>=2)
	      {
		for(unsigned i = 1 ; i < nsam_site ; ++i)
		  denom += 1.0/double(i);
		W += double(nstates-1)/denom;
	      }
	    else if (rep->_totMuts==false &&nstates>=2)
	      {
		for(unsigned i = 1 ; i < nsam_site ; ++i)
		  denom += 1.0/double(i);
		W += 1.0/denom;
	      }
	  }
      }
    return (W);
  }

  double
  PolySNP::ThetaH (void)
  /*!
    Calculate Theta ( = 4Nu) from site homozygosity, a la Fay and Wu (2000).
    This statistic is problematic in general to calculate when there are multiple hits.
    The test requires that the ancestral state (inferred from the outgroup) still be segregating
    in the ingroup.  If that is not true, the site is skipped.\n

    If there are >= 2 derived states inferred, a "missing data" approach is taken.\n
    For example:\n
    Outgroup : \n
    A\n
    Ingroup  : \n
    A\n
    A\n
    A\n
    G\n
    G\n
    T\n
    Gets treated as two sites:\n
    A    A\n
    A    A\n
    A    A\n
    G    N\n
    G    N\n
    N    T\n
    \n
    This keeps the expectation of the statistic equal to \f$\theta\f$, 
    and uses the correct number of derived mutations ovserved in the data.
    \note When using PolySNP, an outgroup is required.  When using PolySIM, 
    which is constructed with a SimData *, an outgroup is not required 
    (as the 0,1 coding of the data refers to ancestral and derived, respectively).
  */
  {
    assert ( rep->_preprocessed );
    if (rep->_NumPoly==0)
      return 0.;
    if( !rep->_haveOutgroup)
      return strtod("NAN",NULL);
    double H = 0.0;
    bool anc_is_present = 0;   //is ancestral state present in the ingroup?

    for (unsigned i = 0;  i < rep->_nsites; ++i)
      {			//iterate over sites
	if ( rep->_derivedCounts[i].second.gap == 0)
	  {
	    unsigned samplesize = rep->_totsam;	//sample size per site
	    unsigned sumDerCounts = 0;
	    sumDerCounts += rep->_derivedCounts[i].second.a;
	    sumDerCounts += rep->_derivedCounts[i].second.g;
	    sumDerCounts += rep->_derivedCounts[i].second.c;
	    sumDerCounts += rep->_derivedCounts[i].second.t;
	    sumDerCounts += rep->_derivedCounts[i].second.zero;
	    sumDerCounts += rep->_derivedCounts[i].second.one;
	    unsigned ancestralCounts = samplesize-sumDerCounts-rep->_derivedCounts[i].second.n;
	    
	    //if the ancestral state is not missing data, and
	    //the sum of the derived counts + and missing data in the ingroup
	    //does not equal the sample size, then the ancestral state is present
	    //at least once in the ingroup, and anc_is_present = true
	    anc_is_present = (rep->_derivedCounts[i].first == true &&
			      ancestralCounts > 0 ) ? true : false; 
	    if (anc_is_present) //ancestral state must be present
	      {
		//number of derived states seen at this site
		int numDer = rep->_derivedCounts[i].second.nStates();	
		if (numDer == 1)
		  {	//simple if there is only one derived state inferred
		    samplesize -= rep->_derivedCounts[i].second.n;
		    double denom = (double(samplesize) *
				    (double(samplesize) - 1.0));
		    H += (rep->_derivedCounts[i].second.a > 0 ) ? 
		      (2.0 * pow (double (rep->_derivedCounts[i].second.a), 2.0)) / denom : 0. ;
		    H += (rep->_derivedCounts[i].second.g > 0 ) ? 
		      (2.0 * pow (double (rep->_derivedCounts[i].second.g), 2.0)) / denom : 0. ;
		    H += (rep->_derivedCounts[i].second.c > 0 ) ? 
		      (2.0 * pow (double (rep->_derivedCounts[i].second.c), 2.0)) / denom : 0. ;
		    H += (rep->_derivedCounts[i].second.t > 0 ) ? 
		      (2.0 * pow (double (rep->_derivedCounts[i].second.t), 2.0)) / denom : 0. ;
		    H += (rep->_derivedCounts[i].second.zero > 0) ? 
		      (2.0 * pow (double (rep->_derivedCounts[i].second.zero), 2.0)) / denom : 0. ;
		    H += (rep->_derivedCounts[i].second.one > 0 ) ? 
		      (2.0 * pow (double (rep->_derivedCounts[i].second.one), 2.0)) / denom : 0. ;
		  }
		else if (numDer == 2 && rep->_haveOutgroup)	//MUST have outgroup--else can't proceed
		  {	//use a "missing data" scheme if there is >1 derived state
		    //iterate over derived states
		    int config[2];
		    unsigned k = 0;
		    if(rep->_derivedCounts[i].second.a > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.a;
		      }
		    if(rep->_derivedCounts[i].second.g > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.g;
		      }
		    if(rep->_derivedCounts[i].second.c > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.c;
		      }
		    if(rep->_derivedCounts[i].second.t > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.t;
		      }
		    if(rep->_derivedCounts[i].second.zero > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.zero;
		      }
		    if(rep->_derivedCounts[i].second.one > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.one;
		      }
		    for (int i = 0; i < 2; ++i)
		      {
			double sample_size_adjust = (i==0) ? double(config[1]) : double(config[0]);
			H += (2.0 *
			      pow (double (config[i]),
				   2.0)) / ((double(samplesize) -
					     sample_size_adjust)
					    * (double(samplesize) -
					       sample_size_adjust
					       - 1.0));
		      }
		  }
	      }
	  }
      }
    return H;
  }

  double
  PolySNP::ThetaL (void)
  /*!
    Calculate Theta ( = 4Nu) from site homozygosity, corresponding to equation 1 in 
    Thornton and Andolfatto (Genetics) "Approximate Bayesian Inference reveals evidence 
    for a recent, severe, bottleneck in a Netherlands population of Drosophila melanogaster," 
    (although we labelled in \f$theta_\eta\f$ in that paper)
    The test requires that the ancestral state (inferred from the outgroup) still be segregating
    in the ingroup.  If that is not true, the site is skipped.\n
	   
    If there are >= 2 derived states inferred, a "missing data" approach is taken.\n
    For example:\n
    Outgroup : \n
    A\n
    Ingroup  : \n
    A\n
    A\n
    A\n
    G\n
    G\n
    T\n
    Gets treated as two sites:\n
    A    A\n
    A    A\n
    A    A\n
    G    N\n
    G    N\n
    N    T\n
    \n
    This keeps the expectation of the statistic equal to \f$\theta\f$, 
    and uses the correct number of derived mutations ovserved in the data.
    \note For sequence data, an outgroup is required.  This requirement is checked by assert()
    @author Joshua Shapiro
  */
  {
    assert(rep->_haveOutgroup==true);
    assert ( rep->_preprocessed );
    double thetal = 0.0;
    if(rep->_NumPoly==0) return thetal;
    bool anc_is_present = 0;   //is ancestral state present in the ingroup?
		  
    for (unsigned i = 0;  i < rep->_nsites; ++i)
      {			//iterate over sites
	if (rep->_derivedCounts[i].first == true && rep->_derivedCounts[i].second.gap == 0)
	  {
	    unsigned samplesize = rep->_totsam;	//sample size per site
	    unsigned sumDerCounts = 0;
	    sumDerCounts += rep->_derivedCounts[i].second.a;
	    sumDerCounts += rep->_derivedCounts[i].second.g;
	    sumDerCounts += rep->_derivedCounts[i].second.c;
	    sumDerCounts += rep->_derivedCounts[i].second.t;
	    sumDerCounts += rep->_derivedCounts[i].second.zero;
	    sumDerCounts += rep->_derivedCounts[i].second.one;
	    unsigned ancestralCounts = samplesize-sumDerCounts-rep->_derivedCounts[i].second.n;
			  
	    //if the ancestral state is not missing data, and
	    //the sum of the derived counts + and missing data in the ingroup
	    //does not equal the sample size, then the ancestral state is present
	    //at least once in the ingroup, and anc_is_present = true
	    anc_is_present = (rep->_derivedCounts[i].first == true &&
			      ancestralCounts > 0 ) ? true : false; 
	    if (anc_is_present) //ancestral state must be present
	      {
		//number of derived states seen at this site
		int numDer = rep->_derivedCounts[i].second.nStates();	
		if (numDer == 1)
		  {	//simple if there is only one derived state inferred
		    samplesize -= rep->_derivedCounts[i].second.n;
		    double denom = (double(samplesize) - 1.0);
		    thetal += (rep->_derivedCounts[i].second.a > 0 ) ? 
		      double (rep->_derivedCounts[i].second.a) / denom : 0. ;
		    thetal += (rep->_derivedCounts[i].second.g > 0 ) ? 
		      double (rep->_derivedCounts[i].second.g) / denom : 0. ;
		    thetal += (rep->_derivedCounts[i].second.c > 0 ) ? 
		      double (rep->_derivedCounts[i].second.c) / denom : 0. ;
		    thetal += (rep->_derivedCounts[i].second.t > 0 ) ? 
		      double (rep->_derivedCounts[i].second.t) / denom : 0. ;
		    thetal += (rep->_derivedCounts[i].second.zero > 0) ? 
		      double (rep->_derivedCounts[i].second.zero) / denom : 0. ;
		    thetal += (rep->_derivedCounts[i].second.one > 0 ) ? 
		      double (rep->_derivedCounts[i].second.one) / denom : 0. ;
		  }
		else if (numDer == 2 && rep->_haveOutgroup)	//MUST have outgroup--else can't proceed
		  {	//use a "missing data" scheme if there is >1 derived state
		    //iterate over derived states
		    int config[2];
		    unsigned k = 0;
		    if(rep->_derivedCounts[i].second.a > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.a;
		      }
		    if(rep->_derivedCounts[i].second.g > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.g;
		      }
		    if(rep->_derivedCounts[i].second.c > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.c;
		      }
		    if(rep->_derivedCounts[i].second.t > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.t;
		      }
		    if(rep->_derivedCounts[i].second.zero > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.zero;
		      }
		    if(rep->_derivedCounts[i].second.one > 0)
		      {
			config[k++] = rep->_derivedCounts[i].second.one;
		      }
		    for (int i = 0; i < 2; ++i)
		      {
			double sample_size_adjust = (i==0) ? double(config[1]) : double(config[0]);
			thetal += (double (config[i])) / (double(samplesize) -
							  sample_size_adjust- 1.0);
		      }
		  }
	      }
	  }
      }
    return thetal;
  }

  unsigned
  PolySNP::NumPoly (void)
  /*!
    \return the number of polymorphic (segregating) sites in data
  */
  {
    assert ( rep->_preprocessed );
    unsigned npoly = 0;
    for (unsigned i = 0;  i < rep->_nsites; ++i)
      {			//iterate over sites
        if (rep->_counts[i].nStates() > 1 && 
	    rep->_counts[i].gap == 0)
          ++npoly;
      }
    return npoly;
  }

  unsigned
  PolySNP::NumMutations (void)
  /*!
    \return the total number of mutations in the data. The number of 
    mutations per site = number of states per site - 1
  */
  {
    assert( rep->_preprocessed );
    unsigned nmut = 0;

    for (unsigned i = 0;  i < rep->_nsites; ++i)
      {			//iterate over sites
        int nstates = (rep->_counts[i].gap==0) ? rep->_counts[i].nStates() : 0;
        if (nstates > 1)
          nmut += nstates - 1;
      }
    return nmut;
  }

  unsigned
  PolySNP::NumSingletons (void)
  /*!
    \return number of polymorphisms that appear once in the data, without respect to ancestral/derived
  */
  {
    assert ( rep->_preprocessed );
    unsigned nsing = 0,nstates;
    for (unsigned i = 0;  i < rep->_nsites; ++i)
      {			//iterate over sites
	unsigned curr_nsing=0,nsam=0;
	nstates = rep->_counts[i].nStates();
	if (rep->_counts[i].gap == 0 && nstates>1)
	  {
	    nsam = rep->_totsam - rep->_counts[i].n;
	    if(nsam==2 && nstates ==2) //if n = 2 and there are 2 states, there must be 1 singleton
	      curr_nsing=1;
	    else
	      {
		curr_nsing += (rep->_counts[i].a == 1) ? 1 : 0;
		curr_nsing += (rep->_counts[i].g == 1) ? 1 : 0;
		curr_nsing += (rep->_counts[i].c == 1) ? 1 : 0;
		curr_nsing += (rep->_counts[i].t == 1) ? 1 : 0;
		curr_nsing += (rep->_counts[i].zero == 1) ? 1 : 0;
		curr_nsing += (rep->_counts[i].one == 1) ? 1 : 0;
	      }
	  }
	nsing += curr_nsing;
      }
    return nsing;
  }


  unsigned
  PolySNP::NumExternalMutations (void)
  /*!
    \return the number of derived singletons.
    \note For sequence data, an outgroup is required. Will return SEQMAXUNSIGNED if that is not the case.  
  */
  {
    if(!rep->_haveOutgroup) return SEQMAXUNSIGNED;
    assert ( rep->_preprocessed );
    int next = 0;
    for (unsigned i = 0;  i < rep->_nsites ; ++i)
      {			//iterate over sites
	unsigned nsam=rep->_totsam;
	unsigned curr_next=0;
	if(rep->_derivedCounts[i].first == true && 
	   rep->_derivedCounts[i].second.gap == 0)
	  {
	    nsam -= rep->_derivedCounts[i].second.n;
	    curr_next += (rep->_derivedCounts[i].second.a == 1) ? 1 : 0;
	    curr_next += (rep->_derivedCounts[i].second.g == 1) ? 1 : 0;
	    curr_next += (rep->_derivedCounts[i].second.c == 1) ? 1 : 0;
	    curr_next += (rep->_derivedCounts[i].second.t == 1) ? 1 : 0;
	    curr_next += (rep->_derivedCounts[i].second.zero == 1) ? 1 : 0;
	    curr_next += (rep->_derivedCounts[i].second.one == 1) ? 1 : 0;
	  }
	next += (nsam>1) ? curr_next : 0;
      }
    return next;
  }


  double
  PolySNP::TajimasD (void)
  /*!
    A common summary of the site frequency spectrum.  
    Proportional to \f$\widehat\theta_\pi-\widehat\theta_W\f$.
    This routine does calculate the denominator of the test statistic.
    \warning statistic undefined if there are untyped SNPs
    \return Tajima's D, or nan if there are no polymorphic sites
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0) return strtod("NAN",NULL);
    double D = 0.0;
    double Pi = ThetaPi ();
    double W = ThetaW ();
    if (fabs(Pi-0.) <= DBL_EPSILON && fabs(W-0.) <= DBL_EPSILON)
      D = 0.0;
    else
      D = (Pi - W) / Dnominator ();
    return D;
  }

  double PolySNP::Hprime (bool likeThorntonAndolfatto)
  /*!
    \return ThetaPi-ThetaH/(~Var(ThetaPi-ThetaH)).  This corresponds to Equation 5 in 
    Thornton and Andolfatto (Genetics) "Approximate Bayesian Inference reveals evidence 
    for a recent, severe, bottleneck in a Netherlands population of Drosophila melanogaster" and Equation 13 of Zeng et al. (2006) Genetics 1431-1439
    \param likeThorntonAndolfatto The calculation of H' requires calculation of 
    \f$\theta^2\f$.  In Thornton and Andolfatto, we simply used \f$\widehat\theta_W^2\f$, 
    which is slightly biased.  By default, this function calculates 
    \f$\theta^2=\frac{S(S-1)}{a_n^2+b_n}\f$, unless this bool is set to false, 
    in which case  \f$\widehat\theta_W^2\f$ is used.
    \note returns nan if there are 0 polymorphic sites
    @author Joshua Shapiro
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0) return strtod("NAN",NULL);
    assert(rep->_haveOutgroup==true);
    double Hpr = 0.0;
    double a = a_sub_n ();
    double b = b_sub_n ();
    double pi = ThetaPi ();
    double theta = ThetaW();
		 
    double thetal = ThetaL();
    double b_n_plus1 = b_sub_n_plus1();
    double S = (rep->_totMuts) ? NumMutations() : NumPoly();
    double thetasq = (likeThorntonAndolfatto == false ) ? S * (S-1)/(a*a + b) : theta*theta;
		  
    double vThetal = 
      (rep->_totsam * theta)/(2.0 * (rep->_totsam - 1.0)) 
      + (2.0 * pow(rep->_totsam/(rep->_totsam - 1.0), 2.0) * (b_n_plus1 - 1.0) - 1.0) * thetasq;
		  
    double vPi = 
      (3.0 * rep->_totsam *(rep->_totsam + 1.0) * theta 
       + 2.0 * ( rep->_totsam * rep->_totsam + rep->_totsam + 3.0) * thetasq  )
      / (9 * rep->_totsam * (rep->_totsam -1.0));
		  
    double cov = 
      ((rep->_totsam + 1.0) / (3.0 * (rep->_totsam - 1.0))) * theta
      +(( 7.0 * rep->_totsam * rep->_totsam + 3.0 * rep->_totsam - 2.0 - 4.0 * rep->_totsam *( rep->_totsam + 1.0) * b_n_plus1)
	/ (2.0 * pow ((rep->_totsam - 1.0), 2.0)))
      * thetasq	 ;
		  
    //    Hpr = pi - omega;
    Hpr = pi - thetal;
    Hpr /= pow ( (vThetal + vPi - 2.0 * cov), 0.5);
    return (Hpr); 
  }


  double
  PolySNP::Dnominator (void)
  /*!
    \warning statistic undefined if there are untyped SNPs
    \return Denominator of Tajima's D, or nan if there are no polymorphic sites
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0) return strtod("NAN",NULL);
    double S = 0.0;
    if (rep->_totMuts)
      {
        S = double (NumMutations ());
      }
    else if (!(rep->_totMuts))
      {
        S = double (NumPoly ());
      }
    double a1, a2, b1, b2, c1, c2, e1, e2;

    a1 = a_sub_n ();
    a2 = b_sub_n ();
    b1 = (rep->_totsam + 1.0) / (3.0 * (rep->_totsam - 1.0));
    b2 = (2.0 * (pow (rep->_totsam, 2.0) 
		 + rep->_totsam + 3.0)) / (9.0 * rep->_totsam *
					   (rep->_totsam - 1.0));
    c1 = b1 - 1.0 / a1;
    c2 = b2 - (rep->_totsam + 2.0) / (a1 * rep->_totsam) + a2 / pow (a1, 2.0);
    e1 = c1 / a1;
    e2 = c2 / (pow (a1, 2.0) + a2);
    double denominator = pow ((e1 * S + e2 * S * (S - 1.0)), 0.5);
    return (denominator);
  }

  void
  PolySNP::DepaulisVeuilleStatistics (void)
  /*!
    Calculate the number of haplotypes in the sample, and haplotype diversity.  
    Unlike Depaulis and Veuille's
    original paper, this routine uses an unbiased calculation of 
    haplotype diversity (i.e. divide by n choose 2).
    \n
    To check if two sequences are unique, Sequence::Comparisons::Different 
    is used, which does not
    allow missing data to result in 2 sequences being considered 
    different (as they would be if you
    simply used the std::string comparison operators == or !=)
  */
  {
    assert ( rep->_preprocessed );
    if (!(rep->_CalculatedDandV))
      {
	if(rep->_NumPoly == 0)
	  {
	    rep->_DVK = 1;
	    rep->_DVH = 0.;
	    return;
	  }
	if (rep->_data->size() > 0)
	  {
	    //step 1 : determine which sequences are unique in the data,
	    //exluding missing data
	    std::set<string,uniqueSeq> unique_haplotypes;
	    if (rep->_haveOutgroup)
	      {
		unique_haplotypes.insert(rep->_data->begin(),
					 rep->_data->begin()+rep->_outgroup);
		unique_haplotypes.insert(rep->_data->begin()+rep->_outgroup+1,
					 rep->_data->end());
	      }
	    else
	      {
		unique_haplotypes.insert(rep->_data->begin(),
					 rep->_data->end());
	      }
	    //now do the real work
	    std::set<string,uniqueSeq>::const_iterator beg = unique_haplotypes.begin(),
	      end = unique_haplotypes.end();
	    rep->_DVK = unsigned(unique_haplotypes.size());
	    while(beg != end)
	      {
		std::ptrdiff_t _count = 0;
		if (rep->_haveOutgroup)
		  {
		    _count += std::count_if(rep->_data->begin(),
					    rep->_data->begin()+rep->_outgroup,
					    boost::bind(notDifferent<string>,
							_1,*beg,false,true));
		    _count += std::count_if(rep->_data->begin()+rep->_outgroup+1,
					    rep->_data->end(),
					    boost::bind(notDifferent<string>,
							_1,*beg,false,true));
		  }
		else
		  {
		    _count += std::count_if(rep->_data->begin(),
					    rep->_data->end(),
					    boost::bind(notDifferent<string>,
							_1,*beg,false,true));
		  }
		rep->_DVH -= pow (double (_count) / rep->_totsam, 2.0);
		++beg;
	      }
	    rep->_DVH *= rep->_totsam / (rep->_totsam - 1.0);
	    rep->_CalculatedDandV = 1;
	  }
      }
  }

  double PolySNP::WallsB(void)
  /*!
    \return Wall's B Statistic. Wall, J. (1999) Genetical Research
    74, pp 65-79
    @author Kevin Thornton
  */
  {
    assert ( rep->_preprocessed );
    if (rep->_calculated_wall_stats == false)
      {
        WallStats();
      }
    return rep->_walls_B;
  }

  void PolySNP::WallStats(void)
  {
    assert ( rep->_preprocessed );
    unsigned S = 0;
    //explicity count # of bi-allelic sites,
    //since that's the proper denominator
    for( std::vector<stateCounter>::const_iterator itr = rep->_counts.begin() ;
	 itr < rep->_counts.end();
	 ++itr)
      {
	if ( itr->nStates() == 2 && itr->gap == 0)
	  ++S;
      }
    if (S > 1)
      {
	std::ptrdiff_t nhap_curr, nhap_left;

	nhap_left = std::ptrdiff_t(SEQMAXUNSIGNED);

	unsigned A = 0;//number of partitions with D' = 1 (see Wall 1999)
	//iterate over sites (actually, adjacent pairs of sites)
	for (unsigned site1 = 0 ; site1 < rep->_nsites-1 ; ++site1)
	  {
	    for(unsigned site2=site1+1 ; site2 < rep->_nsites ; ++site2)
	      {
		if ( rep->_counts[site1].nStates() == 2
		     && rep->_counts[site2].nStates() == 2 )
		  {
		    std::string config;
		    config.resize(2);
		    std::set<string,uniqueSeq> unique_haplotypes;
		    nhap_curr = 0;
		    for (unsigned i = 0 ; i < rep->_nsam ; ++i)
		      {
			if ( (!rep->_haveOutgroup) || (rep->_haveOutgroup && i != rep->_outgroup) )
			  {
			    config[0] = (*rep->_data)[i][site1];
			    config[1] = (*rep->_data)[i][site2];
			    unique_haplotypes.insert(config);
			  }
		      }
		    nhap_curr = unique_haplotypes.size();
		    if(site1==0)
		      {
			if (nhap_curr == 2)
			  {
			    ++rep->_walls_Bprime;
			    ++A;
			  }
		      }
		    else
		      {
			if (nhap_curr == 2)
			  ++rep->_walls_Bprime;
			if (nhap_curr == 2 && nhap_left != 2)
			  ++A;
		      }
		    nhap_left = nhap_curr;
		    site1=site2;
		  }
	      }
	  }
	rep->_walls_B = double(rep->_walls_Bprime)/(double(S-1));
	rep->_walls_Q = (double(rep->_walls_Bprime) + double(A))/(double(S));
      }
    else 
      {
	rep->_walls_B = strtod("NAN",NULL);
	rep->_walls_Bprime = 0;
	rep->_walls_Q = strtod("NAN",NULL);
      }
    rep->_calculated_wall_stats=true;
  }


  unsigned PolySNP::WallsBprime(void)
  /*!
    \return Wall's B' Statistic. Wall, J. (1999) Genetical Research
    74, pp 65-79
    @author Kevin Thornton
  */
  {
    assert ( rep->_preprocessed );
    if (rep->_calculated_wall_stats == false)
      {
        WallStats();
      }
    return rep->_walls_Bprime;
  }

  double PolySNP::WallsQ(void)
  /*!
    \return Wall's Q Statistic. Wall, J. (1999) Genetical Research
    74, pp 65-79
    @author Kevin Thornton
  */
  {
    assert ( rep->_preprocessed );
    if (rep->_calculated_wall_stats == false)
      {
        WallStats();
      }
    return rep->_walls_Q;
  }

  double
  PolySNP::VarPi (void)
  /*!
    Total variance of mean pairwise differences. Tajima in Takahata/Clark book, (13).
    \warning statistic undefined if there are untyped SNPs
  */
  {
    double Pi = ThetaPi ();
    double variance = 3.0 * rep->_totsam * (rep->_totsam + 1.0) * Pi +
      2.0 * (pow (rep->_totsam, 2.0) + rep->_totsam + 3.0) * pow (Pi, 2.0);
    variance /= (11.0 * pow (rep->_totsam, 2.0) - 7.0 * rep->_totsam + 6.0);
    return (variance);
  }

  double
  PolySNP::StochasticVarPi (void)
  /*!
    Stochastic variance of mean pairwise differences. Tajima in Takahata/Clark book, (14).
    \warning statistic undefined if there are untyped SNPs
  */
  {
    double Pi = ThetaPi ();
    double variance = (3.0 * pow (rep->_totsam, 2.0) - 3.0 * rep->_totsam + 2.0) * Pi +
      2.0 * rep->_totsam * (rep->_totsam - 1.0) * pow (Pi, 2.0);
    variance /= (11.0 * pow (rep->_totsam, 2.0) - 7.0 * rep->_totsam + 6.0);
    return (variance);
  }

  double
  PolySNP::SamplingVarPi (void)
  /*!
    Component of variance of mean pairwise differences from sampling. 
    Tajima in Takahata/Clark book, (15)
    \warning statistic undefined if there are untyped SNPs
  */
  {
    double Pi = ThetaPi ();
    double variance =
      2.0 * (3.0 * rep->_totsam - 1.0) * Pi + 2.0 * (2.0 * rep->_totsam +
						     3.0) * pow (Pi, 2.0);
    variance /= (11.0 * pow (rep->_totsam, 2.0) - 7.0 * rep->_totsam + 6.0);
    return (variance);
  }


  double
  PolySNP::VarThetaW (void)
  /*!
    \return Variance of Watterson's Theta (ThetaW()).
    \warning statistic undefined if there are untyped SNPs
  */
  {
    double a1 = a_sub_n ();
    double a2 = b_sub_n ();
    double S = (rep->_totMuts) ? NumMutations() : NumPoly();
    double variance = pow (a1, 2.0) * S + a2 * pow (S, 2.0);
    variance /= pow (a1, 2.0) * (pow (a1, 2.0) + a2);
    return (variance);
  }

  //correct
  double
  PolySNP::FuLiD (void)
  /*!
    \return The Fu and Li (1993) D statistic, or nan if there are no polymorphic sites.
    \note For sequence data, an outgroup is required.  This requirement is checked by assert()
    \warning statistic undefined if there are untyped SNPs
  */
  {
    assert ( rep->_preprocessed );
    //    assert(rep->_haveOutgroup == true);
    if(rep->_NumPoly==0 || !rep->_haveOutgroup) return strtod("NAN",NULL);
    double D = 0.0;
    double ExternalMutations =
      double (NumExternalMutations ());
    double NumMut = double (NumMutations ());
    double a = a_sub_n ();
    double b = b_sub_n ();
    double c = c_sub_n ();
    double vD = 1.0 +
      (pow (a, 2.0) / (b + pow (a, 2.0)) *
       (c - (rep->_totsam + 1.0) / (rep->_totsam - 1.0)));
    double uD = a - 1.0 - vD;
    D = NumMut - a * double (ExternalMutations);
    D /= pow ((uD * NumMut + vD * pow (NumMut, 2.0)), 0.5);
    return (D);
  }


  //correct
  double
  PolySNP::FuLiF (void)
  /*!
    \return Fu and Li (1993) F statistic, or nan if there are no polymorphic sites
    \note For sequence data, an outgroup is required, else undefined
    \warning statistic undefined if there are untyped SNPs
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0 || !rep->_haveOutgroup) return strtod("NAN",NULL);
    double F = 0.0;
    double Pi = ThetaPi ();
    double NumMut = double (NumMutations());
    double ExternalMutations =
      double (NumExternalMutations ());
    double a = a_sub_n ();
    double a_n_plus1 = a_sub_n_plus1 ();
    double b = b_sub_n ();
    double c = c_sub_n ();
    double vF = c + 2.0 * (pow (rep->_totsam, 2.0) + rep->_totsam +
                           3.0) / (9.0 * rep->_totsam * (double (rep->_totsam - 1.0)));
    vF -= (2.0 / (rep->_totsam - 1.0));
    vF /= (pow (a, 2.0) + b);

    double uF = 1.0 + (rep->_totsam + 1.0) / (3.0 * (double (rep->_totsam - 1.0)));
    uF -= 4.0 * ((rep->_totsam + 1.0) / (pow (rep->_totsam - 1.0, 2.0))) *
      (a_n_plus1 - 2.0 * rep->_totsam / (rep->_totsam + 1.0));
    uF /= a;
    uF -= vF;

    F = Pi - ExternalMutations;
    F /= pow (uF * NumMut + vF * pow (NumMut, 2.0), 0.5);
    return (F);
  }

  //correct
  double
  PolySNP::FuLiDStar (void)
  /*!
    \warning statistic undefined if there are untyped SNPs
    \return Fu and Li (1993) D*, or nan if there are no polymorphic sites
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0) return strtod("NAN",NULL);
    double DStar = 0.0;
    double Singletons =
      double (NumSingletons ());
    double NumMut = double (NumMutations ());

    double a = a_sub_n ();
    double b = b_sub_n ();
    double d = d_sub_n ();

    double vD = pow (rep->_totsam / (rep->_totsam - 1.0), 2.0) * b;
    vD += pow (a, 2.0) * d;
    vD -= 2.0 * (rep->_totsam * a * (a + 1.0)) /
      (pow (double (rep->_totsam - 1.0), 2.0));
    vD /= (pow (a, 2.0) + b);

    double uD =
      (rep->_totsam / (rep->_totsam - 1.0)) * (a -
					       (rep->_totsam / (rep->_totsam - 1.0))) - vD;

    DStar = (rep->_totsam / (rep->_totsam - 1.0)) * NumMut - a * double (Singletons);
    DStar /= pow (uD * NumMut + vD * pow (NumMut, 2.0), 0.5);
    return (DStar);
  }

  //correct
  double
  PolySNP::FuLiFStar (void)
  /*!
    Fu and Li (1993) F* statistic. Incorporates correction from
    Simonsen et al.  (1995) Genetics 141: 413, eqn A5.
    \warning statistic undefined if there are untyped SNPs
    \return Fu and Li (1993) F* statistic, or nan if there are no polymorphic sites
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0) return strtod("NAN",NULL);
    double FStar = 0.0;
    double Singletons =
      double (NumSingletons ());
    double Pi = ThetaPi ();
    double NumMut = double (NumMutations ());

    double a = a_sub_n ();
    double a_n_plus1 = a_sub_n_plus1 ();
    double b = b_sub_n ();
    //vF is taken from the correction published by
    //Simonsen et al.  (1995) Genetics 141: 413, eqn A5
    double vF = 2.0 * pow (rep->_totsam, 3.0) + 110.0 * pow (rep->_totsam,
							     2.0) -
      255.0 * rep->_totsam + 153.0;
    vF /= (9.0 * pow (rep->_totsam, 2.0) * (rep->_totsam - 1.0));
    vF += (((2.0 * (rep->_totsam - 1.0) * a) / pow (rep->_totsam, 2.0)) -
           (8.0 * b / rep->_totsam));
    vF /= (pow (a, 2.0) + b);

    double uF =
      (4.0 * pow (rep->_totsam, 2.0) + 19.0 * rep->_totsam + 3.0 -
       12.0 * (rep->_totsam + 1.0) * a_n_plus1);
    uF /= (3.0 * rep->_totsam * (rep->_totsam - 1.0));
    uF /= a;
    uF -= vF;
    FStar = Pi - (((rep->_totsam - 1.0) / rep->_totsam)) * double (Singletons);
    FStar /= pow ((uF * NumMut + vF * pow (NumMut, 2.0)), 0.5);
    return (FStar);
  }

  double
  PolySNP::a_sub_n (void)
  /*!
    \f[a_n=\sum_{i=1}^{i=n-1}\frac{1}{i}.\ \f]
    This is the denominator of Watterson's Theta (see PolySNP::ThetaW)
    \warning statistic undefined if there are untyped SNPs
  */
  {
    assert ( rep->_preprocessed );
    int i;
    double a = 0.0;
    for (i = 1; i < int (rep->_totsam); ++i)
      a += 1. / double (i);
    return a;
  }

  double
  PolySNP::a_sub_n_plus1 (void)
  /*!
    \f[a_{n+1}=\sum_{i=1}^{i=n}\frac{1}{i}\ \f]
    \warning statistic undefined if there are untyped SNPs
  */
  {				//used by Fu and Li tests
    assert ( rep->_preprocessed );
    int i;
    double a = 0.0;
    for (i = 1; i < int (rep->_totsam) + 1; ++i)
      {
        a += 1. / double (i);
      }
    return (a);
  }

  double
  PolySNP::b_sub_n (void)
  /*!
    \f[b_n=\sum_{i=1}^{i=n-1}\frac{1}{i^2}\ \f]
    \warning statistic undefined if there are untyped SNPs
  */
  {				// sum of 1/i^2
    assert ( rep->_preprocessed );
    int i;
    double b = 0.0;
    for (i = 1; i < int (rep->_totsam); ++i)
      b += 1. / (pow (double (i), 2.0));
    return b;
  }

  double
  PolySNP::b_sub_n_plus1(void)
  /*!
    \f[b_n=\sum_{i=1}^{i=n}\frac{1}{i^2}\ \f]
    \warning statistic undefined if there are untyped SNPs
    @author Joshua Shapiro
  */
  {				// sum of 1/i^2
    assert ( rep->_preprocessed );
    int i;
    double b = 0.0;
    for (i = 1; i < int (rep->_totsam) + 1; ++i)
      b += 1. / (pow (double (i), 2.0));
    return b;
  }

  double
  PolySNP::c_sub_n (void)
  /*!
    \f[
    c_n=\left\{\begin{array}{cl}
    1 , & when\ n = 2 \\
    \frac{2 \times (n \times a_n - 2 \times (n-1))}{(n-1) \times (n-2)}, & when\ n > 2 \\
    \end{array}\right.\ 
    \f]
    \warning statistic undefined if there are untyped SNPs
  */
  {				//from Fu and Li 93
    assert ( rep->_preprocessed );
    double c = 0.0, a = a_sub_n ();
    if (fabs(rep->_totsam-2.) <= DBL_EPSILON)
      {
        c = 1.0;
      }
    else
      {
        c = 2.0 * (rep->_totsam * a - 2.0 * (rep->_totsam - 1.0));
        c /= ((rep->_totsam - 1.0) * (rep->_totsam - 2.0));
      }
    return c;
  }

  double
  PolySNP::d_sub_n (void)
  /*!
    \f[\ d_n=\frac{2}{n-1} \times (1.5 - \frac{2 \times a_{n+1}}{n-2} - \frac{1}{n})\ \f]
    \warning statistic undefined if there are untyped SNPs
  */
  {				//from Fu and Li 93
    assert ( rep->_preprocessed );
    double a_n_plus1, c, d;
    a_n_plus1 = a_sub_n_plus1 ();
    c = c_sub_n ();
    d = c + (rep->_totsam - 2.0) / (pow (rep->_totsam - 1.0, 2.0));
    d += (2.0 / (rep->_totsam - 1.0)) *
      (1.5 -  ((2.0 * a_n_plus1 - 3.0) / (rep->_totsam - 2.0)) -  1.0 / rep->_totsam);
    return d;
  }

  double
  PolySNP::DandVH (void)
  /*!
    To check if two sequences are unique, Sequence::Comparisons::Different
    is used, which does not
    allow missing data to result in 2 sequences being 
    considered different (as they would be if you
    simply used thestd::string comparison operators == or !=)
    \return the haplotype diversity of the data.
  */
  {
    if (!(rep->_CalculatedDandV))
      DepaulisVeuilleStatistics ();

    return rep->_DVH;
  }

  unsigned
  PolySNP::DandVK (void)
  /*!
    To check if two sequences are unique, Sequence::Comparisons::Different
    is used, which does not
    allow missing data to result in 2 sequences being considered 
    different (as they would be if you
    simply used the std::string comparison operators == or !=)
    \return number of haplotypes in the sample
  */
  {
    if (!(rep->_CalculatedDandV))
      DepaulisVeuilleStatistics ();

    return rep->_DVK;
  }

  double
  PolySNP::HudsonsC (void)
  /*!
    \return Hudson's (1987) estimator of \f$\rho=4Nc\f$, 
    an estimator of the population recombination rate that 
    depends on the variance of the site frequencies.
    The calculation is made by a call to Recombination::HudsonsC
    \note Will return nan if there are no polymorphic sites
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly==0) return strtod("NAN",NULL);
    return(Recombination::HudsonsC (rep->_data, rep->_haveOutgroup, rep->_outgroup));
  }


  unsigned
  PolySNP::Minrec (void)
  /*!
    \return The minimum number of recombination events observed
    in the sample (Hudson and Kaplan 1985). Will return SEQMAXUNSIGNED 
    if there are < 2 segregating sites. 
    \note Code is a modification of that provided by Jeff Wall
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly<2) return SEQMAXUNSIGNED;
    unsigned a,b,e,numgametes,Rmin=0,x=0;
    bool flag=false;
    
    char c11,c12,c21,c22;
    unsigned states1=0,states2=0;
    
    for (a=x+1 ; a < rep->_nsites ; ++a)
      {
	c11 = c12 = 'Z'; //Z is a dummy value
	//count # states in site a
	states1 = rep->_counts[a].nStates();

	c11 = (c11 == 'Z' && rep->_counts[a].a > 0 ) ? 'A' : 'Z';
	c11 = (c11 == 'Z' && rep->_counts[a].g > 0 ) ? 'G' : c11;
	c11 = (c11 == 'Z' && rep->_counts[a].c > 0 ) ? 'C' : c11;
	c11 = (c11 == 'Z' && rep->_counts[a].t > 0 ) ? 'T' : c11;
	c11 = (c11 == 'Z' && rep->_counts[a].zero > 0 ) ? '0' : c11;
	c11 = (c11 == 'Z' && rep->_counts[a].one > 0 ) ? '1' : c11;
      
	c12 = (c12 == 'Z' && c11 != 'A' && rep->_counts[a].a > 0) ? 'A' : 'Z';
	c12 = (c12 == 'Z' && c11 != 'G' && rep->_counts[a].g > 0) ? 'G' : c12;
	c12 = (c12 == 'Z' && c11 != 'C' && rep->_counts[a].c > 0) ? 'C' : c12;
	c12 = (c12 == 'Z' && c11 != 'T' && rep->_counts[a].t > 0) ? 'T' : c12;
	c12 = (c12 == 'Z' && c11 != '0' && rep->_counts[a].zero > 0) ? '0' : c12;
	c12 = (c12 == 'Z' && c11 != '1' && rep->_counts[a].one > 0) ? '1' : c12;

	for (b = (flag == false) ? x : a-1 ; b < a; ++b)
	  {
	    flag = false;
	    numgametes = 0;
	    c21=c22='Z';
	    states2 = rep->_counts[b].nStates();
	    //need to skip sites with > 2 states
	    if(states1==2&&states2==2)
	      {
		c21 = (c21 == 'Z' && rep->_counts[b].a > 0 ) ? 'A' : 'Z';
		c21 = (c21 == 'Z' && rep->_counts[b].g > 0 ) ? 'G' : c21;
		c21 = (c21 == 'Z' && rep->_counts[b].c > 0 ) ? 'C' : c21;
		c21 = (c21 == 'Z' && rep->_counts[b].t > 0 ) ? 'T' : c21;
		c21 = (c21 == 'Z' && rep->_counts[b].zero > 0 ) ? '0' : c21;
		c21 = (c21 == 'Z' && rep->_counts[b].one > 0 ) ? '1' : c21;
	      
		c22 = (c22 == 'Z' && c21 != 'A' && rep->_counts[b].a > 0) ? 'A' : 'Z';
		c22 = (c22 == 'Z' && c21 != 'G' && rep->_counts[b].g > 0) ? 'G' : c22;
		c22 = (c22 == 'Z' && c21 != 'C' && rep->_counts[b].c > 0) ? 'C' : c22;
		c22 = (c22 == 'Z' && c21 != 'T' && rep->_counts[b].t > 0) ? 'T' : c22;
		c22 = (c22 == 'Z' && c21 != '0' && rep->_counts[b].zero > 0) ? '0' : c22;
		c22 = (c22 == 'Z' && c21 != '1' && rep->_counts[b].one > 0) ? '1' : c22;

		for (e = 0 ; e < rep->_nsam ; ++e)
		  {
		    if (!rep->_haveOutgroup || (rep->_haveOutgroup && e != rep->_outgroup) )
		      if (toupper( (*rep->_data)[e][a] ) == c11  &&
			  toupper( (*rep->_data)[e][b] ) == c21 )
			{
			  ++numgametes;
			  break;
			}
		  }
		for (e = 0 ; e < rep->_nsam ; ++e)
		  {
		    if (!rep->_haveOutgroup || (rep->_haveOutgroup && e != rep->_outgroup) )
		      if (toupper( (*rep->_data)[e][a] ) == c11  &&
			  toupper( (*rep->_data)[e][b] ) == c22 )
			{    
			  ++numgametes;
			  break;
			}
		  }
		for (e = 0 ; e < rep->_nsam ; ++e)
		  {
		    if (!rep->_haveOutgroup || (rep->_haveOutgroup && e != rep->_outgroup) )
		      if (toupper( (*rep->_data)[e][a] ) == c12  &&
			  toupper( (*rep->_data)[e][b] ) == c21 )
			{
			  ++numgametes;
			  break;
			}
		  }
		for (e = 0 ; e < rep->_nsam ; ++e)
		  {
		    if (!rep->_haveOutgroup || (rep->_haveOutgroup && e != rep->_outgroup) )
		      if (toupper( (*rep->_data)[e][a] ) == c12  &&
			  toupper( (*rep->_data)[e][b] ) == c22 )
			{
			  ++numgametes;
			  break;
			}
		  }
		if (numgametes == 4)
		  {
		    ++Rmin;
		    flag = true;
		    break;
		  }
	      }
	  }
	if (flag == true)
	  x=a;
      }
    return Rmin;
  }

  std::vector < std::vector < double > >
  PolySNP::Disequilibrium ( const unsigned & mincount,
			    const double & max_marker_distance)
  /*!
    \return A vector of statistics related to LD and distance in the sample. An empty vector is returned if there are < 2 polymorphic sites in the sample.
    See the documentation for Recombination::Disequilibrium for a 
    description of the return vector.
    \param mincount a frequency filter.  A polymorphism must be present at least \a mincount times in the data
    \note For D and D', the 11 gamete is defined as follows: If no
    outgroup is present, it refers to the genotype of minor alleles at both sites.
    If there is an outgroup, it is based on the genotype of derived alleles at both
    sites.
  */
  {
    assert ( rep->_preprocessed );
    if(rep->_NumPoly<2) return   std::vector < std::vector < double > >();
    return Recombination::Disequilibrium (rep->_data, rep->_haveOutgroup, rep->_outgroup,
					  mincount,max_marker_distance);
  }
}
