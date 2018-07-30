#ifndef SEQUENCE_SUMMSTATS_NSL_HPP__
#define SEQUENCE_SUMMSTATS_NSL_HPP__

#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct nSLiHS
    /// Stores the results of nSL and iHS calculations.
    /// See nsl for details.
    /// \ingroup popgenanalysis
    {
        /// The nSL statistic
        double nsl;
        /// The iHS statistic
        double ihs;
        /// Count of non-reference,
        /// non-missing allele.
        std::int32_t core_count;
    };

    nSLiHS nsl(const VariantMatrix& m, const std::size_t core,
               const std::int8_t refstate);
} // namespace Sequence

#endif
