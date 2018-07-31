/// \file Sequence/summstats/generic.hpp
/// \brief Generic utilities for calculating summary statistics
#ifndef SEQUENCE_SUMMSTATS_GENERIC_HPP
#define SEQUENCE_SUMMSTATS_GENERIC_HPP

#include <unordered_map>
#include <cstdint>

namespace Sequence
{
    /*! \brief Calculate heterozygosity/diversity from count data
     * \param counts a vector counts. 
     * \param nsam the sample size
     * \return diversity = 1 - homozygosity
     */
    double diversity_from_counts(
        const std::unordered_map<std::int32_t, std::int32_t>& counts,
        const std::size_t nsam);
} // namespace Sequence

#endif
