#include <cstdint>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    std::vector<double>
    lhaf(const VariantMatrix &m, const std::int8_t refstate, const double l)
    {
        std::vector<long int> dcounts;
        dcounts.reserve(m.nsites());
        const auto find_nonref = [refstate](const std::int8_t x) {
            return x != refstate && !(x < 0);
        };
        for (std::size_t i = 0; i < m.nsites(); ++i)
            {
                auto r = get_ConstRowView(m, i);
                dcounts.push_back(
                    std::count_if(r.begin(), r.end(), find_nonref));
            }

        // Get the values for each element in the data
        std::vector<double> rv;
        rv.reserve(m.nsam());
        for (std::size_t i = 0; i < m.nsam(); ++i)
            {
                auto c = get_ConstColView(m, i);
                auto j = std::find_if(c.cbegin(), c.cend(), find_nonref);
                double score = 0.0;
                while (j != c.cend())
                    {
                        size_t d2 = static_cast<std::size_t>(
                            std::distance(c.cbegin(), j));
                        score += std::pow(static_cast<double>(dcounts[d2]), l);
                        j = std::find_if(j + 1, c.cend(), find_nonref);
                    }
                rv.push_back(score);
            }
        return rv;
    }
} // namespace Sequence
