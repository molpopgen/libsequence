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
  }
}

#endif
