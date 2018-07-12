#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <cstdint>
#include <unordered_map>

namespace Sequence
{
    std::uint32_t
    nvariable_sites(const VariantMatrix& m)
    {
        std::uint32_t nv = 0;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                std::unordered_map<std::int8_t, std::uint32_t> counts;
                unsigned nmissing = 0;
                for (auto i : site_view)
                    {
                        if (i == VariantMatrix::mask)
                            {
                                throw std::runtime_error(
                                    "invalid state value encountered");
                            }
                        else if (i < 0)
                            {
                                ++nmissing;
                            }
                        else
                            {
                                counts[i]++;
                            }
                    }
                if (counts.size() > 1)
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
                std::unordered_map<std::int8_t, std::uint32_t> counts;
                unsigned nmissing = 0;
                for (auto i : site_view)
                    {
                        if (i == VariantMatrix::mask)
                            {
                                throw std::runtime_error(
                                    "invalid state value encountered");
                            }
                        else if (i < 0)
                            {
                                ++nmissing;
                            }
                        else
                            {
                                counts[i]++;
                            }
                    }
                if (counts.size() > 1)
                    {
                        nv += static_cast<decltype(nv)>(counts.size()) - 1;
                    }
            }
        return nv;
    }
} // namespace Sequence
