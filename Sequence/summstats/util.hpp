/*!
 * \file Sequence/summstats/util.hpp
 * \brief Helper functions for implementing summary statistics
 * 
 * This file must be included directly.  No other header file
 * includes it.
*/
#ifndef SEQUENCE_SUMMSTATS_UTIl_HPP__
#define SEQUENCE_SUMMSTATS_UTIl_HPP__

#include <cstdint>
#include <algorithm>

namespace Sequence
{
    template <typename T>
    inline bool
    all_missing(const T& t)
    /// Returns true if all elements in t encode missing data.
    /// T should be a model of a VariantMatrix, RowView, or ColumnView
    {
        return std::all_of(
            t.begin(), t.end(),
            [](const typename T::value_type v) { return v < 0; });
    }
} // namespace Sequence

#endif
