#include <Sequence/variant_matrix/filtering.hpp>
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>

namespace Sequence
{
    template <typename T>
    std::int32_t
    filter_common_resizable(
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
                                m.position(i)
                                    = std::numeric_limits<double>::quiet_NaN();
                            }
                    }
            }
        if (rv)
            {
                m.resize_capsules();
            }
        return rv;
    }

    template <typename T>
    std::int32_t
    filter_common(
        VariantMatrix &m, const std::function<bool(const T &)> &f,
        const std::function<T(VariantMatrix &, const std::size_t)> &viewmaker,
        std::size_t &dim)
    {
        constexpr auto remove_pos = std::is_same<T, ConstRowView>::value;
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
                                          std::back_inserter(newdata));
                                if (remove_pos)
                                    {
                                        newpos.push_back(m.position(i));
                                    }
                            }
                        else
                            {
                                ++removed;
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
        return filter_common_resizable<RowView>(
            m, f,
            [](VariantMatrix &m, const std::size_t i) {
                return get_RowView(m, i);
            },
            m.nsites());
    }

    std::int32_t
    filter_sites(VariantMatrix &m,
                 const std::function<bool(const ConstRowView &)> &f)
    {
        return filter_common<ConstRowView>(
            m, f,
            [](VariantMatrix &m, const std::size_t i) {
                return get_ConstRowView(m, i);
            },
            m.nsites());
    }

    std::int32_t
    filter_haplotypes(VariantMatrix &m,
                      const std::function<bool(const ColView &)> &f)
    {
        return filter_common_resizable<ColView>(
            m, f,
            [](VariantMatrix &m, const std::size_t i) {
					return get_ColView(m, i);
            },
            m.nsam());
    }

    std::int32_t
    filter_haplotypes(VariantMatrix &m,
                      const std::function<bool(const ConstColView &)> &f)
    {
        return filter_common<ConstColView>(
            m, f,
            [](VariantMatrix &m, const std::size_t i) {
                return get_ConstColView(m, i);
            },
            m.nsam());
    }
} // namespace Sequence
