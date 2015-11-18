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

namespace Sequence
{
  struct variableSiteData
  {
    //! Site positions
    double position;
    //! State counts
    stateCounter counts;
    //! Derived state counts
    stateCounter dcounts;
    //! For each site, is the ancestral state present?
    bool ancState;
    //! Constructo. Copy vs. move left up to caller
    variableSiteData( double, bool, stateCounter,stateCounter );
  };

  template<typename T>
  class variantCounts
  {
  public:
    using data_t = std::vector<variableSiteData>;
  private:
    data_t preprocess(const T & pt, bool haveAnc, unsigned anc, char gapchar)
    {
      data_t rv;
      for( auto i = pt.sbegin() ; i != pt.send() ; ++i )
	{
	  if( haveAnc )
	    {
	      //This is the "counts"
	      auto sc = std::for_each(std::begin(i->second),std::begin(i->second)+anc,stateCounter(gapchar));
	      sc = std::for_each(std::begin(i->second)+anc+1,std::end(i->second),sc);
	      
	      //Now, let's get some 'dcounts'
	      stateCounter sc2(gapchar);
	      auto ancState = *(std::begin(i->second)+anc);
	      bool ancPresent = false;
	      std::for_each(std::begin(i->second),std::begin(i->second)+anc,
			    [&sc2,ancState,&ancPresent](const char & c) {
			      if(c!=ancState) { sc2(c); }
			      else ancPresent=true;
			    });
	      std::for_each(std::begin(i->second)+anc+1,std::end(i->second),
			    [&sc2,ancState,&ancPresent](const char & c) {
			      if(c!=ancState) { sc2(c); }
			      else ancPresent=true;
			    });
	      rv.emplace_back(i->first,ancPresent,std::move(sc),std::move(sc2));
	    }
	  else
	    {
	      //pass empty stateCounter as last arg as anc/derived is taken care of in 'sc'
	      rv.emplace_back( variableSiteData(i->first,false,
						  std::for_each(std::begin(i->second),std::end(i->second),stateCounter()),
						  stateCounter()) );
	    }
	}
      return rv;
    }
  public:
    data_t data;
    using type = T;
    variantCounts(const T & t,bool haveAnc = false, unsigned anc = 0, char gapchar = '-') : data( preprocess(t,haveAnc,anc,gapchar) )
    {
      static_assert(std::is_base_of<PolyTable,T>::value,"T must be derived from Sequence::PolyTable");
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
