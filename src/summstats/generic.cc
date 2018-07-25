#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace Sequence
{
    double
    diversity_from_counts(const std::vector<std::int32_t>& counts,
                          const std::size_t nsam)
    {
        if (counts.empty() || !nsam)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto nmissing

            = std::count_if(counts.begin(), counts.end(),
                            [](decltype(counts[0]) x) { return x < 0; });
        auto nsam_adjusted = nsam - static_cast<decltype(nsam)>(nmissing);
        std::unordered_map<int, int> label_counts;
        for (auto l : counts)
            {
                if (l >= 0)
                    {
                        label_counts[l]++;
                    }
            }
        double hom = 0.0;
        for (auto&& c : label_counts)
            {
                hom += static_cast<double>(c.second * (c.second - 1));
            }
        hom /= static_cast<double>(nsam_adjusted * (nsam_adjusted - 1));
        return 1.0 - hom;
    }
} // namespace Sequence
