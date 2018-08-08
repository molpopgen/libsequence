/// \file Sequence/summstats/thetaw.hpp
/// \brief Watterson's theta
#ifndef SEQUENCE_SUMMSTATS_THETAW_HPP__
#define SEQUENCE_SUMMSTATS_THETAW_HPP__

#include <Sequence/AlleleCountMatrix.hpp>

namespace Sequence
{
    /*! \brief Watterson's theta
     * \param m An AlleleCountMatrix
     * \returns Watterson's theta, a double
     * 
     * \note For a site with \f$k\f$ states,
     * \f$k-1\f$ is added to the number of inferred mutations.
     * In other words, the calculation is based on the total
     * number of mutations.
     *
     * See \cite Watterson1975-ej for details.
     * \ingroup popgenanalysis
     */
    double thetaw(const AlleleCountMatrix& ac);
} // namespace Sequence

#endif
