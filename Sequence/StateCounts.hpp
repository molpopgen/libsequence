#ifndef SEQUENCE_VARIANTMATRIX_STATECOUNTS_HPP__
#define SEQUENCE_VARIANTMATRIX_STATECOUNTS_HPP__

#include "VariantMatrix.hpp"
#include "VariantMatrixViews.hpp"
#include <limits>
#include <vector>

namespace Sequence
{
    struct StateCounts
    /// \brief Track character state occurrence at a site in a VariantMatrix.
    ///
    /// This class keeps track of how many times each character state occurs
    /// at a variable site in a VariantMatrix.  All missing data (negative
    /// state values) are considered equivalent and collapsed into the single
    /// missing value of -1.
    ///
    /// When constructed, the sample size at a site is considered to be the
    /// sum of the number of occurrences of all non-missing states.
    ///
    /// \ingroup variantmatrix
    {
        static constexpr VariantMatrix::value_type max_allele
            = std::numeric_limits<VariantMatrix::value_type>::max();
        /// Keep track of (state, count) pairs
        std::vector<std::int32_t> counts;
        /// The max allelic value seen
        std::size_t max_allele_idx;
        /// The sample size at this site.  Excluded missing data.
        std::uint32_t n;
        /// The reference state for this site.  Needed for certain summary
        /// statistics. Default is -1 (missing).
        std::int8_t refstate;

        /// Construct with a ConstRowView and a reference state, which defaults
        /// to 0.
        StateCounts(const std::int8_t refstate_);
        StateCounts();
        void operator()(ConstRowView &);
        void operator()(const RowView &);
    };

    /// Create a vector of StateCounts from a VariantMatrix.
    /// If `refstates` is not empty and differs in length
    /// from `m.nsites`, then `std::invalid_argument` is thrown.
    /// \ingroup variantmatrix
    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m,
                           const std::vector<std::int8_t>& refstates);
    /// Create a vector of StateCounts with a specific reference state
    /// used for all sites
    /// \ingroup variantmatrix
    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m, const std::int8_t refstate);
    /// Create a vector of StateCounts with a reference state of -1
    /// used for all sites
    /// \ingroup variantmatrix
    std::vector<StateCounts> process_variable_sites(const VariantMatrix& m);
} // namespace Sequence

#endif
