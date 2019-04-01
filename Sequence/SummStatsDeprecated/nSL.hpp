#ifndef __SEQUENCE_SUMMSTATS_NSL_HPP__
#define __SEQUENCE_SUMMSTATS_NSL_HPP__

/* \file nSL.hpp
   @brief The nSL statistic of doi: 10.1093/molbev/msu077
*/

#include <Sequence/PolySIM.hpp>
#include <unordered_map>
#include <tuple>
#include <cstdint>

namespace Sequence
{
    /*!
      The nSL statistic of Ferrer-Admetlla et al. doi: 10.1093/molbev/msu077.
      \param core The index of the "focal/core" SNP
      \param d An object of type Sequence::SimData
      \param gmap The positions of every marker in d on the genetic map.  If
      std::unordered_map<double,double>() is passed,
      iHS is calculated using SNP positions.
      \return nSL and iHs, with the latter as defined in doi:
      10.1093/molbev/msu077.
      \note This routine was validated by comparing to code provided by
      Ferrer-Admetlla et al.
      \warning The use of 'gmap' is untested.
      \ingroup popgenanalysis
     */
    std::pair<double, double>
    nSL(const std::size_t &core, const SimData &d,
        const std::unordered_map<double, double> &gmap
        = std::unordered_map<double, double>())__attribute__ ((deprecated));
}
#endif
