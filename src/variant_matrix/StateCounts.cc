#include <Sequence/StateCounts.hpp>
#include <algorithm>
#include <stdexcept>

inline static decltype(Sequence::StateCounts::n)
fill_counts(decltype(Sequence::StateCounts::counts)& counts,
            const Sequence::ConstRowView& r)
{

    decltype(Sequence::StateCounts::n) n{ 0u };
    for (auto i : r)
        {
            if (i == Sequence::VariantMatrix::mask)
                {
                    throw std::invalid_argument("reserved value encountered");
                }
            if (i < 0)
                i = -1;
            auto itr = counts.find(i);
            if (itr == counts.end())
                {
                    counts.insert(std::make_pair(i, 1));
                }
            else
                {
                    itr->second++;
                }
        }
    for (const auto& ci : counts)
        {
            if (ci.first >= 0)
                {
                    n += ci.second;
                }
        }
    return n;
}

namespace Sequence
{
    StateCounts::StateCounts(const ConstRowView& r)
        : counts{}, n{ 0u }, refstate(-1)
    {
        n = fill_counts(counts, r);
    }

    StateCounts::StateCounts(const ConstRowView& r,
                             const std::int8_t refstate_)
        : counts{}, n{ 0u }, refstate(refstate_)
    {
        n = fill_counts(counts, r);
    }

    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m,
                           const std::vector<std::int8_t>& refstates)
    {
        if (std::any_of(
                refstates.begin(), refstates.end(),
                [](const std::int8_t x) { return x == VariantMatrix::mask; }))
            {
                throw std::invalid_argument("reserved value encountered");
            }
        if (refstates.size() != m.nsites)
            {
                throw std::invalid_argument("refstates.size() != m.nsites");
            }
        std::vector<StateCounts> rv;
        rv.reserve(m.nsites);
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                rv.emplace_back(get_ConstRowView(m, i), refstates[i]);
            }
        return rv;
    }

    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m, const std::int8_t refstate)
    {
        if (refstate == Sequence::VariantMatrix::mask)
            {
                throw std::invalid_argument("reserved value encountered");
            }
        std::vector<StateCounts> rv;
        rv.reserve(m.nsites);
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                rv.emplace_back(get_ConstRowView(m, i), refstate);
            }
        return rv;
    }

    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m)
    {
        return process_variable_sites(m, -1);
    }
} // namespace Sequence
