#include <cstdint>
#include <vector>
#include <algorithm>
#include <Sequence/summstats/ld.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    TwoLocusCounts::TwoLocusCounts(std::int8_t i_, std::int8_t j_, int n_)
        : i{ i_ }, j{ j_ }, n{ n_ }
    {
    }

    std::vector<TwoLocusCounts>
    two_locus_haplotype_counts(const VariantMatrix& m, std::size_t sitei,
                               const std::size_t sitej,
                               const bool skip_missing)
    {
        auto ri = get_ConstRowView(m, sitei);
        auto rj = get_ConstRowView(m, sitej);
        std::vector<TwoLocusCounts> rv;
        for (auto i = ri.begin(), j = rj.begin(); i < ri.end(); ++i, ++j)
            {
                if (!skip_missing || ((*i < 0 || *j < 0) && !skip_missing))
                    {
                        auto exists
                            = std::find_if(rv.begin(), rv.end(),
                                           [i, j](const TwoLocusCounts& t) {
                                               return t.i == *i && t.j == *j;
                                           });
                        if (exists == rv.end())
                            {
                                rv.emplace_back(*i, *j, 1);
                            }
                        else
                            {
                                exists->n++;
                            }
                    }
            }
        return rv;
    }
} // namespace Sequence
