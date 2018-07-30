#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace Sequence
{
    double
    diversity_from_counts(
        const std::unordered_map<std::int32_t, std::int32_t>& counts,
        const std::size_t nsam)
    {
        if (counts.empty() || !nsam)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        double hom = 0.0;
        for (auto&& c : counts)
            {
                hom += static_cast<double>(c.second * (c.second - 1));
            }
        hom /= static_cast<double>(nsam * (nsam - 1));
        return 1.0 - hom;
    }
} // namespace Sequence
