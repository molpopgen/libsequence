#include <cstdint>
#include <algorithm>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    std::uint32_t
    nvariable_sites(const VariantMatrix& m)
    {
        std::uint32_t nv = 0;
        StateCounts counts;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                counts(site_view);
                auto nstates = std::count_if(
                    counts.counts.begin(), counts.counts.end(),
                    [](const std::int32_t x) { return x > 0; });
                if (nstates > 1)
                    {
                        ++nv;
                    }
            }
        return nv;
    }

    std::uint32_t
    nbiallelic_sites(const VariantMatrix& m)
    {
        std::uint32_t nv = 0;
        StateCounts counts;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                counts(site_view);

                auto nstates = std::count_if(
                    counts.counts.begin(), counts.counts.end(),
                    [](const std::int32_t x) { return x > 0; });
                if (nstates == 2)
                    {
                        ++nv;
                    }
            }
        return nv;
    }

    std::uint32_t
    total_number_of_mutations(const VariantMatrix& m)
    {
        std::uint32_t nv = 0;
        StateCounts counts;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                counts(site_view);

                auto nstates = std::count_if(
                    counts.counts.begin(), counts.counts.end(),
                    [](const std::int32_t x) { return x > 0; });
                if (nstates > 1)
                    {
                        nv += static_cast<decltype(nv)>(nstates) - 1;
                    }
            }
        return nv;
    }
} // namespace Sequence
