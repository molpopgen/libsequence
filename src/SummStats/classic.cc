#include <iterator>
#include <iostream>
#include <algorithm>
#include <set>
#include <cmath>
#include <Sequence/SummStats/classic_details.hpp>
#include <Sequence/Comparisons.hpp>

namespace
{
  struct uniqueSeq : public std::binary_function<std::string,std::string,bool>
 {
    inline bool operator()(const std::string & l, const std::string  & r) const
    {
      //use Sequence::Different to prevent missing sites 
      //causing 2 sequences to be labelled as distinct
      //return (  Different(l,r) && std::lexicographical_compare(l.begin(),l.end(),r.begin(),r.end(),lt_nocase()) );
      return (  Sequence::Different(l,r) && std::lexicographical_compare(l.begin(),l.end(),r.begin(),r.end(),
									 [](const char & __a, const char __b)
									 {
									   return (std::toupper(static_cast<unsigned char>(__a)) 
										   < std::toupper(static_cast<unsigned char>(__b)));
									 }) );
    }
  };
}

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

    unsigned singletons_details( const std::vector<variableSiteData> & c, unsigned, std::true_type )
    {
      unsigned rv = 0;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++ i )
	{
	  i += (i->counts.one==1||i->counts.zero==1) ? 1 : 0;
	}
      return rv;
    }

    unsigned dsingletons_details( const std::vector<variableSiteData> & c, unsigned, std::true_type )
    {
      unsigned rv = 0;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++ i )
	{
	  i += (i->counts.one==1) ? 1 : 0;
	}
      return rv;
    }
    
    unsigned singletons_details( const std::vector<variableSiteData> & c, unsigned nsam, std::false_type )
    {
      unsigned rv = 0;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++ i )
	{
	  if(!i->counts.gap && (nsam-i->counts.n > 1))
	    {
	      rv += (i->counts.a == 1) ? 1u : 0u;
	      rv += (i->counts.g == 1) ? 1u : 0u;
	      rv += (i->counts.c == 1) ? 1u : 0u;
	      rv += (i->counts.t == 1) ? 1u : 0u;
	      rv += (i->counts.zero == 1) ? 1u : 0u;
	      rv += (i->counts.one == 1) ? 1u : 0u;
	    }
	}
      return rv;
    }
    unsigned dsingletons_details( const std::vector<variableSiteData> & c, unsigned nsam, std::false_type )
    {
      unsigned rv = 0;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++ i )
	{
	  if(!i->counts.gap&&(nsam-i->counts.n > 1))
	    {
	      rv += (i->dcounts.a == 1) ? 1u : 0u;
	      rv += (i->dcounts.g == 1) ? 1u : 0u;
	      rv += (i->dcounts.c == 1) ? 1u : 0u;
	      rv += (i->dcounts.t == 1) ? 1u : 0u;
	      rv += (i->dcounts.zero == 1) ? 1u : 0u;
	      rv += (i->dcounts.one == 1) ? 1u : 0u;
	    }
	}
      return rv;
    }
  }

  std::pair<unsigned,double> hapstats(const PolyTable & t, const bool haveAnc, const unsigned anc)
  {
    if(t.empty()) return std::make_pair(1u,0.);
    std::set<std::string,uniqueSeq> unique_haplotypes;
    if(haveAnc)
      {
	unique_haplotypes.insert(std::begin(t),std::begin(t)+anc);
	unique_haplotypes.insert(std::begin(t)+anc+1,std::end(t));
      }
    else
      {
	unique_haplotypes.insert(std::begin(t),std::end(t));
      }
    double hdiv=1.0;
    auto beg = std::begin(unique_haplotypes);
    auto end = std::end(unique_haplotypes);
    unsigned nsam = unsigned(t.size() - unsigned(haveAnc));
    for( auto i =beg ; i != end ; ++i )
      {
	double c = 0.;
	if(haveAnc)
	  {
	    c += double(std::count_if(std::begin(t),std::begin(t)+anc,[i](const std::string & __s) {
		  return !Different(__s,*i,false,true);
		}));
	    c += double(std::count_if(std::begin(t)+anc+1,std::end(t),[i](const std::string & __s) {
		  return !Different(__s,*i,false,true);
		}));
	  }
	else
	  {
	    c += double(std::count_if(std::begin(t),std::end(t),[i](const std::string & __s) {
		  return !Different(__s,*i,false,true);
		}));	    
	  }
	hdiv -= std::pow(c/nsam,2.0);
      }
    hdiv *= nsam/(nsam-1.);
    return std::make_pair(unique_haplotypes.size(),hdiv);
  }
  
  std::pair<unsigned,double> hapstats(const SimData & t)
  {
    if(t.empty()) return std::make_pair(1u,0.);
    std::set<std::string> s(t.begin(),t.end());
    double hdiv = 1.0;
    double nsam=double(t.size());
    for(auto i = std::begin(s) ; i != std::end(s) ; ++i )
      {
	double c = double(std::count(std::begin(s),std::end(s),*i));
	hdiv -= std::pow(c/nsam,2.0);
      }
    hdiv *= nsam/(nsam-1.);
    return std::make_pair(s.size(),hdiv);
  }
  
}
