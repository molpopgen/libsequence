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

#include <Sequence/HKA.hpp>
#include <cmath>
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

using std::isfinite;

namespace 
{
  double Cn(unsigned n)
  {
    double cn=0.;
    for(unsigned i=1;i<n;++i)
      {
	cn += 1./double(i);
      }
    return cn;
  }

  double Cnsq(unsigned n)
  {
    double cn=0.;
    for(unsigned i=1;i<n;++i)
      {
	cn += 1./double(i*i);
      }
    return cn;
  }
}

namespace Sequence
{
  HKAdata::HKAdata() : SA(0),SB(0),D(0),nA(0),nB(0)
		  /*!
		    default constructor
		  */
  {
  }

  HKAdata::HKAdata(const HKAdata & d) : 
    SA(d.SA),SB(d.SB),D(d.D),nA(d.nA),nB(d.nB)
		  /*!
		    copy constructor
		  */
  {
  }
  HKAdata::HKAdata(unsigned sa,unsigned sb,double d,
	  unsigned na,unsigned nb) :
    SA(sa),SB(sb),D(d),nA(na),nB(nb)
		  /*!
		    \param sa Num. polymorphic sites in species a
		    \param sb Num. polymorphic sites in species b
		    \param d  Divergence between species a and b (per locus)
		    \param na sample size for species a
		    \param nb sample size for species b
		  */

  {
  }

  HKAresults::HKAresults (const std::vector<double> & _thetas,
			  const std::vector< chisq_tuple > & _chisquareds,
			  const double & _fhat, const double & _That,
			  const double & _xsq,
			  const double & _xsqA,
			  const double & _xsqB)
    : thetas(_thetas),chisquareds(_chisquareds),fhat(_fhat),
      That(_That),xsq(_xsq),xsqA(_xsqA),xsqB(_xsqB)
    /*!
      \param _thetas a vector containing the estimats of theta for each locus
      \param _chisquareds a vector of the per-locus contributions to the test statistic
      \param _fhat an estimate of the relative Ne of species b to species a
      \param _That an estimate of the divergence time, on the scale of Ne in species a.
      \param _xsq the Chi-squared statistic for the HKA test
      \param _xsqA  The HKA chi-squared statistic, only considering species A    
      \param _xsqB  The HKA chi-squared statistic, only considering species B    
    */
  {
  }

  HKAresults calcHKA ( const std::vector< HKAdata > & data )
  /*!
    Performs the calculations necessary for the HKA test.
    \param data a vector of HKAdata objects
    \return an object of type multiLocusHKAparams.  
    \note The thetas vector in the return object contains theta values in the same order as the loci
    in \a data.
    \ingroup popgen
  */
  {
    double sumD=0.;
    //get some summary stats that we need
    double sumTheta = 0.;
    double fhat = 0.;
    double cna,cnb=0.;
    double that=0.;
    unsigned n=0,sumSb=0;
    for(unsigned i=0;i<data.size();++i)
      {
	sumD += data[i].D;   // sum of divergence accross loci
	cna = Cn(data[i].nA);
	cnb += (data[i].nB > 1) ? Cn(data[i].nB) : 0.;
	sumTheta += double(data[i].SA)/cna;
	sumSb += data[i].SB;
	n += (data[i].nB>1)?1:0;
	//double f = (data[i].nB>1 && data[i].SA>0) ? (double(data[i].SB)*cna)/(double(data[i].SA)*cnb): 1.;
	//	n += (data[i].nB>=1 && data[i].SA>0) ? 1 : 0;
	//	fhat += f;
      }
    cnb /= n;
    fhat = ( isfinite(cnb) && cnb > 0. ) ? double(sumSb)/(sumTheta*cnb) : 1.;
    that = double(sumD)/sumTheta - (0.5*(1.+fhat));
    std::vector<double> thetas(data.size());
    std::vector< HKAresults::chisq_tuple > chisquareds;
    double xsq = 0.,xsqA = 0., xsqB = 0;
    //now, get the theta estimates for each locus,
    //and count up the xsq statistic

    double ESA,VSA,ESB,VSB,ED,VD;
    for(unsigned i=0;i<data.size();++i)
      {
	cna = Cn(data[i].nA);
	cnb = Cn(data[i].nB);
	thetas[i] = double(double(data[i].SA) +
			   double(data[i].SB)+ 
			   data[i].D)/
	  (that+0.5*(1+fhat)+cna+fhat*cnb);
	
	ESA = thetas[i]*cna; 
	VSA = ESA + thetas[i]*thetas[i]*Cnsq(data[i].nA);
	ESB =  fhat*thetas[i]*cnb;
	VSB =  ESB + fhat*fhat*thetas[i]*thetas[i]*Cnsq(data[i].nB);
	ED = thetas[i]*(that + 0.5*(1.+fhat));
	VD = ED + (thetas[i]*0.5*(1.+fhat))*(thetas[i]*0.5*(1.+fhat));
	double xsq_t=0.,xsqA_t=0.,xsqB_t=0.,xsqpoly=0.;
	const double spA = (double(data[i].SA)-ESA)*(double(data[i].SA)-ESA)/VSA;
	xsq_t += spA;
	xsqA_t += spA;
	xsqpoly = xsqA_t;
	const double sp2 = (double(data[i].SB)-ESB)*(double(data[i].SB)-ESB)/VSB;
	if (isfinite(sp2)) 
	  {
	    xsq_t += sp2;
	    xsqpoly += sp2;
	  }
	xsqB_t += sp2; //this will be not finite if n=1 in species 2, which is the desired behavior
	const double d = (data[i].D-ED)*(data[i].D-ED)/VD;
	xsq_t += d;
	chisquareds.push_back( boost::make_tuple(xsqpoly,d,xsqA_t,xsqB_t) );
	xsq += xsq_t;
	xsqA_t += d;
	xsqB_t += d;
	xsqA += xsqA_t;
	xsqB += xsqB_t;
      }
    return HKAresults( thetas,chisquareds,fhat,that,xsq,xsqA,xsqB );
  }
}
