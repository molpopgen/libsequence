#ifndef __SEQUENCE_VARIANT_COUNTS__
#define __SEQUENCE_VARIANT_COUNTS__

#include <memory>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <type_traits>
#include <Sequence/stateCounter.hpp>
#include <Sequence/SimData.hpp>
#include <iostream>
#include <limits>

namespace Sequence
{
  struct variableSiteData
  {
    //! Site positions
    double pos;
    //! State counts
    stateCounter counts;
    //! Derived state counts
    stateCounter dcounts;
    //! For each site, is the ancestral state present?
    bool ancState;
    //! Constructor. Copy vs. move left up to caller
    variableSiteData( double, bool, stateCounter &&,stateCounter && );
  };

  template<typename T>
  class variantCounts
  {
  public:
    using data_t = std::vector<variableSiteData>;
  private:
    void preprocess(const T & pt, bool haveAnc, unsigned anc, char gapchar)
    {
      auto nsites = pt.numsites();
      auto nchrom = pt.size();
      data.reserve(nsites);
      if(!haveAnc)
	{
	  for(unsigned i=0;i<nsites;++i)
	    {
	      stateCounter sc(gapchar);
	      for(unsigned j=0;j<nchrom;++j)
		{
		  sc(pt[j][i]);
		}
	      data.emplace_back(pt.position(i),false,std::move(sc),stateCounter(gapchar));
	    }
	}
      else
	{
	  for(unsigned i=0;i<nsites;++i)
	    {
	      bool ancState = false;
	      auto ancstate = pt[anc][i];
	      stateCounter sc(gapchar),scd(gapchar);
	      for(unsigned j=0;j<nchrom;++j)
		{
		  if( j != anc )
		    {
		      auto c_ = pt[j][i];
		      if (c_ != ancstate) scd(c_);
		      else ancState=true;
		      sc(c_);
		    }
		}
	      data.emplace_back(pt.position(i),ancState,std::move(sc),std::move(scd));
	    }
	}
    }
  public:
    data_t data;
    double pi;
    const unsigned nsam;
    unsigned npoly;
    using type = T;
    variantCounts(const T & t,bool haveAnc = false,
		  unsigned anc = 0, char gapchar = '-') :
      data( std::vector<variableSiteData>() ),
      pi(std::numeric_limits<double>::quiet_NaN()),
      nsam(unsigned(t.size()-unsigned(haveAnc))),
      npoly(std::numeric_limits<unsigned>::max())
    {
      static_assert(std::is_base_of<PolyTable,T>::value,"T must be derived from Sequence::PolyTable");
      this->preprocess(t,haveAnc,anc,gapchar);
    }
  };

  template<typename T>
  variantCounts<T> make_variantCounts(const T & t, bool haveAnc = false, unsigned anc = 0, char gapchar = '-' )
  {
    return variantCounts<T>(t,haveAnc,anc,gapchar);
  }

  template<>
  variantCounts<SimData> make_variantCounts(const SimData & sd, bool haveAnc, unsigned anc, char gapchar);
}


#endif
