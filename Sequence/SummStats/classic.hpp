#ifndef SEQUENCE_SUMMSTATS_CLASSIC_HPP
#define SEQUENCE_SUMMSTATS_CLASSIC_HPP
/*!
  "Classic" summaries of nucleotide variability.

  This file makes Sequence::PolySNP/Sequence::PolySIM obsolete.
 */

#include <type_traits>
#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/classic_details.hpp>
#include <Sequence/SummStats/variantCounts.hpp>

namespace Sequence
{
  template<typename T>
  inline double thetapi(const T & t)
  {
    return details::thetapi_details(t.data,t.nsam,typename std::is_same<typename T::type,SimData>::type());
  }
}

#endif
