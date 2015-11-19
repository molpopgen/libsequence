#include <Sequence/SummStats/classic_details.hpp>
#include <iterator>
#include <iostream>
namespace Sequence
{
  namespace details
  {
    double a_sub_n(unsigned n)
    {
      double rv=0.;
      for(unsigned i=0;i<n;++i) rv += 1./double(i);
      return rv;
    }
    
    double thetapi_details( const std::vector<variableSiteData> & c, unsigned nsam_, std::true_type )
    {
      double Pi = 0.0, nsam = double(nsam_);
      for( auto i = std::begin(c) ; i != std::end(c) ; ++i )
	{
	  Pi += (2.0 * i->counts.one *(nsam-i->counts.one)) / (nsam*(nsam-1.));
	}
      return Pi;
    }
    
    double thetapi_details( const std::vector<variableSiteData> & c, unsigned nsam_, std::false_type )
    {
      double Pi = 0.0;
      for(auto i = std::begin(c) ; i != std::end(c) ; ++i )
	{			//iterate over sites
	  if ( i->counts.gap == 0 &&
	       i->counts.nStates() > 1 )
	    {
	      unsigned samplesize = nsam_;
	      samplesize -= i->counts.n; //adjust sample size for missing data
	      if (samplesize > 1)
		{
		  double SSH = 0.0;	//sum of site homozygosity
		  double denom = (double(samplesize)* (double(samplesize) - 1.0));
		  SSH += (i->counts.a > 0) ? double(i->counts.a) * 
		    double (i->counts.a-1) /denom : 0. ;
		  SSH += (i->counts.g > 0) ? double(i->counts.g) * 
		    double (i->counts.g-1) /denom : 0. ;
		  SSH += (i->counts.c > 0) ? double(i->counts.c) * 
		    double (i->counts.c-1) /denom : 0. ;
		  SSH += (i->counts.t > 0) ? double(i->counts.t) * 
		    double (i->counts.t-1) /denom : 0. ;
		  SSH += (i->counts.zero > 0) ? double(i->counts.zero) * 
		    double (i->counts.zero-1) /denom : 0. ;
		  SSH += (i->counts.one > 0) ? double(i->counts.one) * 
		    double (i->counts.one-1) /denom : 0. ;
		  Pi += (1.0 - SSH);
		}
	    }
	}
      return Pi;
    }

    double thetaw_details( const std::vector<variableSiteData> & c, unsigned nsam, bool, std::true_type )
    {
      return double(npoly_details(c))/a_sub_n(nsam);
    }
    
    double thetaw_details( const std::vector<variableSiteData> & c, unsigned nsam, bool totMuts, std::false_type )
    {
      double w = 0.;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++i )
	{
	  if(! i->counts.gap )
	    {
	      auto ns = i->counts.nStates();
	      if(ns>1)
		{
		  w += (totMuts) ? double(ns-1)/a_sub_n(nsam-i->counts.n) : 1.0/a_sub_n(nsam-i->counts.n);
		}
	    }
	}
      return w;
    }
  }
}
