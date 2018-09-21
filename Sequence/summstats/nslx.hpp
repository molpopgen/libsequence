/// \file Sequence/summstats/nslx.hpp
/// \brief nSL and iHS
#ifndef SEQUENCE_SUMMSTATS_NSLX_HPP
#define SEQUENCE_SUMMSTATS_NSLX_HPP

#include <vector>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>
#include "nSLiHS.hpp"

namespace Sequence
{
    /*! \brief A variation on nSL/iHS 
     * \param m A VariantMatrix
     * \param refstate The ancestral state
     * \param x Non-reference allele count
     *
     * \return vector of nSLiHS
     *
     * This variant on nSL only allows suffix lengths
     * to be broken by variants where the derived 
     * (non-refstate) allele is present <= \a x times.
     *
     * When \x is 1, this statistic is a proxy for the 
     * SDS score of \cite Field2016-so.
     */
    std::vector<nSLiHS> nslx(const VariantMatrix& m,
                             const std::int8_t refstate, const int x);
} // namespace Sequence

#endif
