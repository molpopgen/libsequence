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
        auto pb
            = std::lower_bound(m.positions.begin(), m.positions.end(), beg);
        auto pe = std::upper_bound(pb, m.positions.end(), end);
        decltype(m.data) data;
        decltype(m.positions) pos;
        if (pb == m.positions.end())
            {
                return VariantMatrix(std::move(data), std::move(pos));
            }

        pos.assign(pb, pe);
        for (auto i = pb; i < pe; ++i)
            {
                auto v = get_ConstRowView(
                    m, static_cast<std::size_t>(
                           std::distance(m.positions.begin(), i)));
                data.insert(data.end(),v.begin(),v.end());
            }
        return VariantMatrix(std::move(data), std::move(pos));
    }

} // namespace Sequence
