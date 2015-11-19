#ifndef SEQUENCE_SUMMSTATS_CLASSIC_DETAILS_HPP
#define  SEQUENCE_SUMMSTATS_CLASSIC_DETAILS_HPP

#include <type_traits>
#include <Sequence/SummStats/variantCounts.hpp>

namespace Sequence
{
  namespace details
  {
    double thetapi_details( const std::vector<variableSiteData> & c, unsigned nsam, std::true_type );
    double thetapi_details( const std::vector<variableSiteData> & c, unsigned nsam, std::false_type );
    double thetaw_details( const std::vector<variableSiteData> & c, unsigned nsam, bool totMuts,std::true_type );
    double thetaw_details( const std::vector<variableSiteData> & c, unsigned nsam, bool totMuts,std::false_type );
    inline unsigned npoly_details( const std::vector<variableSiteData> & c )
    {
      unsigned np=0;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++i )
      {
	if(i->counts.nStates()>1&&i->counts.gap==0) ++np;
      }
      return np;
    }
    inline unsigned nmuts_details( const std::vector<variableSiteData> & c )
    {
      unsigned rv = 0;
      for( auto i = std::begin(c) ; i != std::end(c) ; ++i )
	{
	  rv += ( !i->counts.gap ) ? (i->counts.nStates()-1) : 0;
	}
      return rv;
    }
  }
}

#endif
