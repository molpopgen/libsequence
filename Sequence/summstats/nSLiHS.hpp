#ifndef SEQUENCE_SUMMSTATS_NSLIHS_HPP
#define SEQUENCE_SUMMSTATS_NSLIHS_HPP

#include <cstdint>

namespace Sequence
{
    struct nSLiHS
    /// Stores the results of nSL and iHS calculations.
    /// See Sequence::nsl for details.
    ///
    /// \note This type is usually forward-declared in other headers,
    /// meaning this header will need inclusion in relevant translation
    /// units.
    ///
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
} // namespace Sequence

#endif
