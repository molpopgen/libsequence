#include <Sequence/StateCounts.hpp>
#include <stdexcept>

namespace Sequence
{
    StateCounts::StateCounts(const ConstRowView& r,
                             const std::int8_t refstate_)
        : counts{}, n{ 0u }, refstate(refstate_)
    {
        for (auto i : r)
            {
                if (i == VariantMatrix::mask)
                    {
                        throw std::invalid_argument(
                            "reserved value encountered");
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
        for(const auto & ci : counts)
        {
            if(ci.first >= 0)
            {
                n += ci.second;
            }
        }
    }

    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m,
                           const std::vector<std::int8_t>& refstates)
    {
        if (!refstates.empty() && refstates.size() != m.nsites)
            {
                throw std::invalid_argument("refstates.size() != m.nsites");
            }
        std::vector<StateCounts> rv;
        rv.reserve(m.nsites);
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                if (refstates.empty())
                    {
                        rv.emplace_back(get_ConstRowView(m, i));
                    }
                else
                    {
                        rv.emplace_back(get_ConstRowView(m, i), refstates[i]);
                    }
            }
        return rv;
    }
}
