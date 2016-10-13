#ifndef __SEQUENCE_SUMMSTATS_NSL_HPP__
#define __SEQUENCE_SUMMSTATS_NSL_HPP__

/* \file nSL.hpp
   @brief The nSL statistic of doi: 10.1093/molbev/msu077
*/

#include <Sequence/PolySIM.hpp>
#include <unordered_map>

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
    nSL(const unsigned &core, const SimData &d,
        const std::unordered_map<double, double> &gmap
        = std::unordered_map<double, double>());

    /*!
      Calculate max. abs value of standardized nSL and iHS, with the latter as
      defined by Ferrer-Admetella et al.
      \param d An object of type Sequence::SimData
      \param minfreq Exclude mutations with minor allele frequency < minfreq.
      \param binsize The size of frequency bins.
      \param gmap The positions of every marker in d on the genetic map. If
      std::unordered_map<double,double>() is passed,
      iHS is calculated using SNP positions.
      \return maximum absolute value of standardized nSL and iHS, with the
      latter as defined by Ferrer-Admetella et al.
      \warning The use of 'gmap' is untested.
      \item The first member of the return value is nSL, the second is iHS
      \ingroup popgenanalysis
    */
    std::pair<double, double>
    snSL(const SimData &d, const double minfreq, const double binsize,
         const std::unordered_map<double, double> &gmap
         = std::unordered_map<double, double>());
}
#endif
