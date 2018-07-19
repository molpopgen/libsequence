/// \file Sequence/summstats/allele_counts.hpp
/// \brief Count alleles at variable sites.
#ifndef SEQUENCE_SUMMSTATS_ALLELE_COUNTS_HPP__
#define SEQUENCE_SUMMSTATS_ALLELE_COUNTS_HPP__

#include <vector>
#include <utility>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct AlleleCounts
    /// Tracks the number of states at a site
    /// \ingroup popgenanalysis
    {
        /// Number of non-missing states
        int nstates;
        /// Number of samples with missing states
        int nmissing;
    };

    /*! \brief Count number of alleles at each site
     * \param m A VariantMatrix
     * \ingroup popgenanalysis
     */
    std::vector<AlleleCounts> allele_counts(const VariantMatrix& m);

    /*! \brief Count number of non-reference alleles at each site
     * \param m A VariantMatrix
     * \param m refstate The reference state for all sites.
     * \ingroup popgenanalysis
     */
    std::vector<AlleleCounts>
    non_reference_allele_counts(const VariantMatrix& m,
                                const std::int8_t refstate);

    /*! \brief Count number of non-reference alleles at each site
     * \param m A VariantMatrix
     * \param m refstate The reference state at each site.
     * \ingroup popgenanalysis
     */
    std::vector<AlleleCounts>
    non_reference_allele_counts(const VariantMatrix& m,
                                const std::vector<std::int8_t>& refstates);
} // namespace Sequence
#endif
