#include <Sequence/NonOwningCapsules.hpp>
#include <Sequence/variant_matrix/windows.hpp>

namespace Sequence
{
    VariantMatrix
    make_window(const VariantMatrix& m, const double beg, const double end)
    {
        if (end < beg)
            {
                throw std::invalid_argument("end must be >= beg");
            }
        auto pb = std::lower_bound(m.pbegin(), m.pend(), beg);
        auto pe = std::upper_bound(pb, m.pend(), end);
        std::vector<int8_t> data;
        std::vector<double> pos;
        if (pb == m.pend())
            {
                return VariantMatrix(std::move(data), std::move(pos));
            }

        pos.assign(pb, pe);
        for (auto i = pb; i < pe; ++i)
            {
                auto v = get_ConstRowView(
                    m, static_cast<std::size_t>(std::distance(m.pbegin(), i)));
                data.insert(data.end(), v.begin(), v.end());
            }
        return VariantMatrix(std::move(data), std::move(pos));
    }

    VariantMatrix
    make_slice(const VariantMatrix& m, const double beg, const double end,
               const std::size_t i, const std::size_t j)
    {
        if (end < beg)
            {
                throw std::invalid_argument("end must be >= beg");
            }
        if (!(j > i))
            {
                throw std::invalid_argument("i must be < j");
            }
        if (i >= m.nsam() || j >= m.nsam())
            {
                throw std::invalid_argument("slice indexes out of range");
            }
        auto pb = std::lower_bound(m.pbegin(), m.pend(), beg);
        auto pe = std::upper_bound(pb, m.pend(), end);
        std::vector<std::int8_t> data;
        std::vector<double> pos;
        if (pb == m.pend() || i == j)
            {
                return VariantMatrix(std::move(data), std::move(pos));
            }
        pos.assign(pb, pe);
        for (auto pi = pb; pi < pe; ++pi)
            {
                auto v = get_ConstRowView(
                    m,
                    static_cast<std::size_t>(std::distance(m.pbegin(), pi)));
                data.insert(data.end(), v.begin() + i, v.begin() + j);
            }
        return VariantMatrix(std::move(data), std::move(pos));
    }
} // namespace Sequence
