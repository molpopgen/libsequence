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
  inline double thetaw(T & t, bool totMuts = true)
  {
    return details::thetaw_details(t.data,t.nsam,totMuts,typename std::is_same<typename T::type,SimData>::type());
  }

  template<typename T>
  inline unsigned npoly(T & t)
  {
    if(t.npoly != std::numeric_limits<unsigned>::max()) return t.npoly;
    t.npoly=details::npoly_details(t.data);
    return t.npoly;
  }

  template<typename T>
  inline unsigned nmuts(T & t)
  {
    return details::nmuts_details(t.data);
  }

  template<typename T>
  inline unsigned singletons(const T & t)
  {
    return details::singletons_details(t.data,typename std::is_same<typename T::type,SimData>::type());
  }
}

#endif
