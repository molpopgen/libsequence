#include <vector>
#include <utility>
#include <cstdint>
#include <Sequence/summstats/allele_counts.hpp>
#include <Sequence/AlleleCountMatrix.hpp>

using count_type = std::vector<Sequence::AlleleCounts>;

namespace
{
    template <typename T>
    typename count_type::value_type
    add_counts(const T& row, const std::size_t nsam, const bool nonref,
               const std::int8_t refstate)
    {
        if (nonref && refstate == -1)
            {
                return count_type::value_type{ -1, -1 };
            }
        std::size_t refindex = static_cast<std::size_t>(refstate);
        if (nonref
            && refindex >= static_cast<std::size_t>(row.second - row.first))
            {
                throw std::invalid_argument("reference state out of range");
            }
        typename count_type::value_type rv{ 0, 0 };
        int nnon_missing = 0;
        for (auto i = row.first; i != row.second; ++i)
            {
                if (*i > 0)
                    {
                        nnon_missing += *i;
                        if (!nonref
                            || (nonref
                                && static_cast<std::size_t>(i - row.first)
                                       != refindex))
                            {
                                ++rv.nstates;
                            }
                    }
            }
        rv.nmissing = static_cast<int>(nsam) - nnon_missing;
        return rv;
    }
} // namespace

namespace Sequence
{
    count_type
    allele_counts(const AlleleCountMatrix& m)
    {
        count_type rv;
        for (std::size_t i = 0; i < m.nrow; ++i)
            {
                auto r = m.row(i);
                rv.emplace_back(add_counts(r, m.nsam, false, -1));
            }
        return rv;
    }

    count_type
    non_reference_allele_counts(const AlleleCountMatrix& m,
                                const std::vector<std::int8_t>& refstates)
    {
        if (refstates.size() != m.nrow)
            {
                throw std::invalid_argument("number of reference states does "
                                            "not equal number of sites");
            }
        count_type rv;
        for (std::size_t i = 0; i < m.nrow; ++i)
            {
                auto r = m.row(i);
                rv.emplace_back(add_counts(r, m.nsam, true, refstates[i]));
            }
        return rv;
    }

    count_type
    non_reference_allele_counts(const AlleleCountMatrix& m,
                                const std::int8_t refstate)
    {
        count_type rv;
        for (std::size_t i = 0; i < m.nrow; ++i)
            {
                auto r = m.row(i);
                rv.emplace_back(add_counts(r, m.nsam, true, refstate));
            }
        return rv;
    }
} // namespace Sequence
