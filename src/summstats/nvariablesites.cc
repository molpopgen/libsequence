#include <cstdint>
#include <algorithm>
#include <Sequence/AlleleCountMatrix.hpp>

namespace Sequence
{
    std::uint32_t
    nvariable_sites(const AlleleCountMatrix& m)
    {
        std::uint32_t nv = 0;
        for (std::size_t site = 0; site < m.nrow; ++site)
            {
                auto r = m.row(site);
                auto nstates
                    = std::count_if(r.first, r.second,
                                    [](const AlleleCountMatrix::value_type c) {
                                        return c > 0;
                                    });
                if (nstates > 1)
                    {
                        ++nv;
                    }
            }
        return nv;
    }

    std::uint32_t
    nbiallelic_sites(const AlleleCountMatrix& m)
    {
        std::uint32_t nv = 0;
        for (std::size_t site = 0; site < m.nrow; ++site)
            {
                auto r = m.row(site);
                auto nstates
                    = std::count_if(r.first, r.second,
                                    [](const AlleleCountMatrix::value_type c) {
                                        return c > 0;
                                    });
                if (nstates == 2)
                    {
                        ++nv;
                    }
            }
        return nv;
    }

    std::uint32_t
    total_number_of_mutations(const AlleleCountMatrix& m)
    {
        std::uint32_t nv = 0;
        for (std::size_t site = 0; site < m.nrow; ++site)
            {
                auto r = m.row(site);
                auto nstates
                    = std::count_if(r.first, r.second,
                                    [](const AlleleCountMatrix::value_type c) {
                                        return c > 0;
                                    });
                if (nstates > 1)
                    {
                        nv += static_cast<decltype(nv)>(nstates) - 1;
                    }
            }
        return nv;
    }
} // namespace Sequence
