#include <cstdint>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>

static inline double
updateH(const Sequence::StateCounts& c, const std::int8_t refstate)
{
    double H = 0.0;
    double nnm1 = static_cast<double>(c.n * (c.n - 1));
    for (auto& x : c.counts)
        {
            if (!(x.first < 0) && x.first != refstate)
                {
                    H += std::pow(x.second, 2.0);
                }
        }
    H /= nnm1;
    return H;
}

namespace Sequence
{
    double
    thetah(const VariantMatrix& m, const std::int8_t refstate)
    {
        if (refstate < 0)
            {
                throw std::invalid_argument(
                    "all reference states are encoded as missing");
            }
        double H = 0.0;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                StateCounts counts(r, refstate);
                H += updateH(counts, refstate);
            }
        return H;
    }

    double
    thetah(const VariantMatrix& m, const std::vector<std::int8_t>& refstates)
    {
        if (std::all_of(refstates.begin(), refstates.end(),
                        [](const std::int8_t x) { return x < 0; }))
            {
                throw std::invalid_argument(
                    "all reference states are encoded as missing");
            }
        double H = 0.0;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                StateCounts counts(r, refstates[i]);
                H += updateH(counts, refstates[i]);
            }
        return H;
    }
} // namespace Sequence
