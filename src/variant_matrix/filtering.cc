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
        std::vector<std::size_t> removed_indexes;
        for (std::size_t i = 0; i < dim; ++i)
            {
                auto view = viewmaker(m, i);
                if (f(view))
                    {
                        ++rv;
                        removed_indexes.push_back(i);
                    }
            }
        if (rv)
            {
                std::vector<double> newpos;
                std::vector<std::int8_t> newdata;
                std::size_t removed = 0;
                for (std::size_t i = 0; i < dim; ++i)
                    {
                        if (removed < removed_indexes.size()
                            && i != removed_indexes[removed])
                            {
                                auto view = viewmaker(m, i);
                                std::copy(view.begin(), view.end(),
                                          end(newdata));
                                if (remove_pos)
                                    {
                                        newpos.push_back(m.positions[i]);
                                    }
                            }
                    }
                VariantMatrix v(std::move(newdata), std::move(newpos));
				swap(m, v);
            }
        return rv;
    }

    std::int32_t
    filter_sites(VariantMatrix &m,
                 const std::function<bool(const RowView &)> &f)
    {
        return filter_common<RowView>(
            m, f,
            [](VariantMatrix &m, const std::size_t i) {
                return get_RowView(m, i);
            },
            m.nsites);
    }

    std::int32_t
    filter_haplotypes(VariantMatrix &m,
                      const std::function<bool(const ColView &)> &f)
    {
        return filter_common<ColView>(
            m, f,
            [](VariantMatrix &m, const std::size_t i) {
                return get_ColView(m, i);
            },
            m.nsam);
    }
} // namespace Sequence
