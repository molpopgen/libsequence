#include <cstdint>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    std::uint32_t
    nvariable_sites(const VariantMatrix& m)
    {
        std::uint32_t nv = 0;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                StateCounts counts(site_view);
                counts.counts.erase(-1);
                if (counts.counts.size() > 1)
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
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                StateCounts counts(site_view);
                counts.counts.erase(-1);
                if (counts.counts.size() > 1)
                    {
                        nv += static_cast<decltype(nv)>(counts.counts.size())
                              - 1;
                    }
            }
        return nv;
    }
} // namespace Sequence
