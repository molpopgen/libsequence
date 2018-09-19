#ifndef SEQUENCE_SUMMSTATS_NSL_COMMON_HPP
#define SEQUENCE_SUMMSTATS_NSL_COMMON_HPP

// These functions are not exported.
// They are used internally.

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <Sequence/summstats/nSLiHS.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    namespace summstats_details
    {
        struct suffix_edges
        {
            std::int64_t left, right;
            suffix_edges() : left(-1), right(-1) {}
        };

        static void
        update_counts(double nsl_values[2], double ihs_values[2],
                      int counts[2], const std::size_t nsites,
                      const std::vector<double>& positions,
                      const std::size_t index, const std::int64_t left,
                      const std::int64_t right)
        {
            if (left >= 0 && static_cast<std::size_t>(right) < nsites)
                //Then there are SNPs differentiating
                //i and j within the region
                {
                    nsl_values[index] += static_cast<double>(right - left);
                    //TODO: check if we need to add one?
                    ihs_values[index]
                        += positions[static_cast<std::size_t>(right)]
                           - positions[static_cast<std::size_t>(left)];
                    counts[index]++;
                }
        }

        inline nSLiHS
        get_stat(const ConstRowView& core_view, const std::int8_t refstate,
                 const double nsl_values[2], const double ihs_values[2],
                 const int counts[2])
        {

            double nSL_den = nsl_values[0] / static_cast<double>(counts[0]);
            double nSL_num = nsl_values[1] / static_cast<double>(counts[1]);
            double iHS_den = ihs_values[0] / static_cast<double>(counts[0]);
            double iHS_num = ihs_values[1] / static_cast<double>(counts[1]);
            auto nonrefcount = static_cast<std::int32_t>(
                std::count_if(core_view.begin(), core_view.end(),
                              [refstate](const std::int8_t i) {
                                  return i >= 0 && i != refstate;
                              }));
            return nSLiHS{ std::log(nSL_num) - std::log(nSL_den),
                           std::log(iHS_num) - std::log(iHS_den),
                           nonrefcount };
        }
    } // namespace summstats_details
} // namespace Sequence

#endif
