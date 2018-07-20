#ifndef SEQUENCE_SUMMSTATS_NSL_HPP__
#define SEQUENCE_SUMMSTATS_NSL_HPP__

#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct nSL
    {
        /// The nSL statistic
        double nsl;
        /// The iHS statistic
        double ihs;
        /// Count of non-reference,
        /// non-missing allele.
        std::int32_t core_count;
    };

    nSL nsl(const VariantMatrix& m, const std::size_t core,
            const std::int8_t refstate);
} // namespace Sequence

#endif
