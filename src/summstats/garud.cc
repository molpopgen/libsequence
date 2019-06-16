#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <Sequence/summstats/garud.hpp>
#include <Sequence/summstats/generic.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/summstats/classics.hpp>

namespace Sequence
{
    GarudStats::GarudStats()
        : H1(1.), H12(std::numeric_limits<double>::quiet_NaN()),
          H2H1(std::numeric_limits<double>::quiet_NaN())
    {
    }

    GarudStats::GarudStats(const double __h1, const double __h12,
                           const double __h2h1)
        : H1(__h1), H12(__h12), H2H1(__h2h1)
    {
    }

    GarudStats
    garud_statistics(const VariantMatrix& m)
    {
        GarudStats rv;
        if (m.empty() || !m.nsam())
            {
                return rv;
            }
        // Although one of the stats is haplotype diversity,
        // w
        auto labels = label_haplotypes(m);
        std::unordered_map<std::int32_t, std::int32_t> counts;
        std::size_t nmissing = 0;
        for (auto l : labels)
            {
                if (l < 0)
                    {
                        ++nmissing;
                    }
                else
                    {
                        counts[l]++;
                    }
            }
        if (counts.size() < 2)
            {
                return rv;
            }
        rv.H1 = 1.0 - diversity_from_counts(counts, m.nsam() - nmissing);
        std::vector<std::pair<std::int32_t, std::int32_t>> vcounts(
            counts.begin(), counts.end());
        std::sort(vcounts.begin(), vcounts.end(),
                  [](const std::pair<std::int32_t, std::int32_t>& a,
                     const std::pair<std::int32_t, std::int32_t>& b) {
                      return a.second > b.second;
                  });
        double nsam
            = static_cast<double>(m.nsam()) - static_cast<double>(nmissing);
        rv.H12 = rv.H1
                 + 2. * static_cast<double>(vcounts[0].second)
                       * static_cast<double>(vcounts[1].second)
                       / (nsam * (nsam - 1.0));
        rv.H2H1 = (rv.H1
                   - static_cast<double>(vcounts[0].second
                                         * (vcounts[0].second - 1))
                         / (nsam * (nsam - 1)))
                  / rv.H1;
        return rv;
    }
} // namespace Sequence
