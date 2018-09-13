#include <Sequence/StateCounts.hpp>
#include <algorithm>
#include <stdexcept>

namespace Sequence
{
    StateCounts::StateCounts()
        : counts(std::vector<std::int32_t>(max_allele + 1, 0)),
          max_allele_idx{ 0 }, n{ 0u }, refstate(-1)
    {
    }

    StateCounts::StateCounts(const std::int8_t refstate_)
        : counts(std::vector<std::int32_t>(max_allele + 1, 0)),
          max_allele_idx{ 0 }, n{ 0u }, refstate(refstate_)
    {
    }

    void
    StateCounts::operator()(ConstRowView& row)
    {
        std::fill(counts.begin(), counts.begin() + max_allele_idx + 1, 0);
        n = 0;
        max_allele_idx = 0;
        for (std::size_t i = 0; i < row.size(); ++i)
            {
                auto ri = row[i];
                if (ri >= 0)
                    {
                        ++n;
                        ++counts[static_cast<std::size_t>(ri)];
                        max_allele_idx = std::max(
                            max_allele_idx, static_cast<std::size_t>(ri));
                    }
                else if (ri == VariantMatrix::mask)
                    {
                        throw std::invalid_argument(
                            "reserved value encountered");
                    }
            }
    }

    void
    StateCounts::operator()(const RowView& row)
    {
        std::fill(counts.begin(), counts.begin() + max_allele_idx + 1, 0);
        n = 0;
        max_allele_idx = 0;
        for (std::size_t i = 0; i < row.size(); ++i)
            {
                auto ri = row[i];
                if (ri >= 0)
                    {
                        ++n;
                        ++counts[static_cast<std::size_t>(ri)];
                        max_allele_idx = std::max(
                            max_allele_idx, static_cast<std::size_t>(ri));
                    }
                else if (ri == VariantMatrix::mask)
                    {
                        throw std::invalid_argument(
                            "reserved value encountered");
                    }
            }
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
                StateCounts c(refstates[i]);
                auto r = get_RowView(m, i);
                c(r);
                rv.emplace_back(std::move(c));
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
                StateCounts c(refstate);
                auto r = get_RowView(m, i);
                c(r);
                rv.emplace_back(std::move(c));
            }
        return rv;
    }

    std::vector<StateCounts>
    process_variable_sites(const VariantMatrix& m)
    {
        return process_variable_sites(m, -1);
    }
} // namespace Sequence
