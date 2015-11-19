#ifndef SEQUENCE_SUMMSTATS_CLASSIC_HPP
#define SEQUENCE_SUMMSTATS_CLASSIC_HPP
/*!
  "Classic" summaries of nucleotide variability.

  This file makes Sequence::PolySNP/Sequence::PolySIM obsolete.
 */

#include <cmath>
#include <type_traits>
#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/classic_details.hpp>
#include <Sequence/SummStats/variantCounts.hpp>

namespace Sequence
{
  template<typename T>
  inline double thetapi(T & t)
  {
    if(std::isfinite(t.pi)) return t.pi;
    t.pi = details::thetapi_details(t.data,t.nsam,typename std::is_same<typename T::type,SimData>::type());
    return t.pi;
  }

  template<typename T>
  inline unsigned npoly(T & t)
  {
    if(t.npoly != std::numeric_limits<unsigned>::max()) return t.npoly;
    t.npoly=0;
    for( auto i = std::begin(t.data) ; i != std::end(t.data) ; ++i )
      {
	if(i->counts.nStates()>1&&i->counts.gap==0) ++t.npoly;
      }
    return t.npoly;
  }

  template<typename T>
  inline unsigned nmuts(T & t)
  {
    unsigned rv = 0;
    for( auto i = std::begin(t.data) ; i != std::end(t.data) ; ++i )
      {
	rv += ( !i->counts.gap ) ? (i->counts.nStates()-1) : 0;
      }
    return rv;
  }
}

#endif
