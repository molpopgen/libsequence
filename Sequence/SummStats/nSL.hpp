#ifndef __SEQUENCE_SUMMSTATS_NSL_HPP__
#define __SEQUENCE_SUMMSTATS_NSL_HPP__

/* \file nSL.hpp
   @brief The nSL statistic of doi: 10.1093/molbev/msu077
*/

#include <Sequence/PolySIM.hpp>

namespace Sequence {

  /*
    The nSL statistic of Ferrer-Admetlla et al. doi: 10.1093/molbev/msu077.
    \param core The index of the "focal/core" SNP
    \param d An object of type Sequence::SimData
    \return The nSL statistic with respect to core. If the statistic is 
    not defined for this site, a non-finite value will be returned.
    \note This routine was validated by comparing to code provided by
    Ferrer-Admetlla et al.
   */
  double nSL(const unsigned & core,
	     const SimData & d);
}
#endif
