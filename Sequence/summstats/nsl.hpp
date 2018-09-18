/// \file Sequence/summstats/nsl.hpp
/// \brief nSL and iHS
#ifndef SEQUENCE_SUMMSTATS_NSL_HPP__
#define SEQUENCE_SUMMSTATS_NSL_HPP__

#include <vector>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct nSLiHS
    /// Stores the results of nSL and iHS calculations.
    /// See Sequence::nsl for details.
    /// \ingroup popgenanalysis
    {
        /// The nSL statistic \cite Ferrer-Admetlla2014-wa
        double nsl;
        /// The iHS statistic, calculated according to \cite Ferrer-Admetlla2014-wa
        double ihs;
        /// Count of non-reference,
        /// non-missing allele.
        std::int32_t core_count;
    };

    nSLiHS nsl(const VariantMatrix& m, const std::size_t core,
               const std::int8_t refstate);

    std::vector<nSLiHS> nsl(const VariantMatrix& m,
                            const std::int8_t refstate);
} // namespace Sequence

#endif
