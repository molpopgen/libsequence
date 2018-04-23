#ifndef SEQUENCE_VARIANTMATRIX_STATECOUNTS_HPP__
#define SEQUENCE_VARIANTMATRIX_STATECOUNTS_HPP__

#include "VariantMatrixViews.hpp"
#include <unordered_map>
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
        /// Keep track of (state, count) pairs
        std::unordered_map<std::int8_t, std::uint32_t> counts;
        /// The sample size at this site.  Excluded missing data.
        std::uint32_t n;
        /// The reference state for this site.  Needed for certain summary
        /// statistics. Default is -1 (missing).
        std::int8_t refstate;

		/// Construct with a ConstRowView and a reference state, which defaults
		/// to 0.
        StateCounts(const ConstRowView& r, const std::int8_t refstate_ = 0);
    };

	/// Create a vector of StateCounts from a VariantMatrix.
	/// If you wish to specify a reference state separately for each site,
	/// do so by passing in `refstates`, which much have the same length
	/// as `m.nsites`.  If `refstates` is not empty and differs in length
	/// from `m.nsites`, then `std::invalid_argument` is thrown.
	/// \ingroup variantmatrix
    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m,
                           const std::vector<std::int8_t>& refstates
                           = std::vector<int8_t>());
}

#endif
