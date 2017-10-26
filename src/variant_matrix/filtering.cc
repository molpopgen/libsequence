#include <Sequence/variant_matrix/filtering.hpp>
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>

namespace Sequence
{
    template <typename T>
    std::int32_t
    filter_common(
        VariantMatrix &m, const std::function<bool(const T &)> &f,
        const std::function<T(VariantMatrix &, const std::size_t)> &viewmaker,
        std::size_t &dim)
    {
        constexpr auto remove_pos = std::is_same<T, RowView>::value;
        std::int32_t rv = 0;
        for (std::size_t i = 0; i < dim; ++i)
            {
                auto view = viewmaker(m, i);
                if (f(view))
                    {
                        ++rv;
                        std::transform(
                            view.begin(), view.end(), view.begin(),
                            [](double) { return VariantMatrix::mask; });
                        if (remove_pos)
                            {
                                m.positions[i]
                                    = std::numeric_limits<double>::quiet_NaN();
                            }
                    }
            }
        if (rv)
            {
                if (remove_pos)
                    {
                        m.positions.erase(
                            std::remove_if(
                                m.positions.begin(), m.positions.end(),
                                [](double d) { return std::isnan(d); }),
                            m.positions.end());
                    }
                m.data.erase(std::remove(m.data.begin(), m.data.end(),
                                         VariantMatrix::mask),
                             m.data.end());
                dim -= rv;
            }
        return rv;
    }

    std::int32_t
    filter_sites(VariantMatrix &m,
                 const std::function<bool(const RowView &)> &f)
    {
        return filter_common<RowView>(
            m, f, [](VariantMatrix &m,
                     const std::size_t i) { return get_RowView(m, i); },
            m.nsites);
    }

    std::int32_t
    filter_haplotypes(VariantMatrix &m,
                      const std::function<bool(const ColView &)> &f)
    {
        return filter_common<ColView>(
            m, f, [](VariantMatrix &m,
                     const std::size_t i) { return get_ColView(m, i); },
            m.nsam);
    }
}
