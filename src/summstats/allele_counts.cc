#include <vector>
#include <utility>
#include <cstdint>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

using count_type = std::vector<std::pair<std::int32_t, std::int32_t>>;

namespace
{
    typename count_type::value_type
    add_counts(const Sequence::StateCounts& counts, const bool nonref)
    {
        if (nonref && counts.refstate == -1)
            {
                return count_type::value_type{ -1, -1 };
            }
        typename count_type::value_type rv{ 0, 0 };
        for (auto& c : counts.counts)
            {
                if (c.first < 0)
                    {
                        ++rv.second;
                    }
                else if (!nonref || (nonref && c.first != counts.refstate))
                    {
                        ++rv.first;
                    }
            }
        return rv;
    }
} // namespace

namespace Sequence
{
    count_type
    allele_counts(const VariantMatrix& m)
    {
        count_type rv;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                StateCounts counts(r);
                rv.emplace_back(add_counts(counts, false));
            }
        return rv;
    }

    count_type
    non_reference_allele_counts(const VariantMatrix& m,
                                const std::vector<std::int8_t>& refstates)
    {
        if (refstates.size() != m.nsites)
            {
                throw std::invalid_argument("number of reference states does "
                                            "not equal number of sites");
            }
        count_type rv;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                StateCounts counts(r, refstates[i]);
                rv.emplace_back(add_counts(counts, true));
            }
        return rv;
    }

    count_type
    non_reference_allele_counts(const VariantMatrix& m,
                                const std::int8_t refstate)
    {
        count_type rv;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                StateCounts counts(r, refstate);
                rv.emplace_back(add_counts(counts, true));
            }
        return rv;
    }
} // namespace Sequence
