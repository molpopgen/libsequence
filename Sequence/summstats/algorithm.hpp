#ifndef SEQUENCE_SUMMSTATS_ALGORITHM_HPP
#define SEQUENCE_SUMMSTATS_ALGORITHM_HPP

#include <cstdint>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    namespace sstats_algo
    {
        template <typename F>
        inline void
        aggregate_sites(const VariantMatrix& m, const F& f,
                        const std::int8_t refstate)
        {
            StateCounts c(refstate);
            for (std::size_t i = 0; i < m.nsites; ++i)
                {
                    auto r = get_ConstRowView(m, i);
                    c(r);
                    f(c);
                }
        }
    } // namespace sstats_algo
} // namespace Sequence

#endif
