#ifndef SEQUENCE_SUMMSTATS_GENERIC_HPP
#define SEQUENCE_SUMMSTATS_GENERIC_HPP

#include <vector>
#include <cstdint>

namespace Sequence
{
    double diversity_from_counts(const std::vector<std::int32_t>& counts,
                                 const std::size_t nsam);
}

#endif
