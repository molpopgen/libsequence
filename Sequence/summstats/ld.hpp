#ifndef SEQUENCE_SUMMSTATS_LD_HPP__
#define SEQUENCE_SUMMSTATS_LD_HPP__

#include <cstdint>
#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct TwoLocusCounts
    {
        std::int8_t i, j;
        int n;
    };

    std::vector<TwoLocusCounts>
    two_locus_haplotype_counts(const VariantMatrix& m, std::size_t sitei,
                               const std::size_t sitej,
                               const bool skip_missing);
} // namespace Sequence

#endif
