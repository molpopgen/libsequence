#include <algorithm>
#include <Sequence/summstats/nsl.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <iostream>
#include "nsl_common.hpp"

namespace
{
    inline bool
    update_edge_matrix(const Sequence::VariantMatrix& m,
                       const std::vector<std::int64_t>& xtons,
                       const std::pair<std::int64_t, std::int64_t> flanks,
                       const Sequence::ConstRowView& core_view,
                       const Sequence::ConstColView& hapi,
                       std::vector<std::int64_t>& edges,
                       const std::size_t core, const std::size_t i,
                       const std::size_t j)
    // edge updating for nslx
    // precondition: !xtons.empty()
    {
        bool rv = false;
        std::size_t lindex = j * m.nsam + i;
        if (core_view[i] == core_view[j] && !(core_view[i] < 0))
            {
                auto hapj = get_ConstColView(m, j);
                std::int64_t left_edge = edges[lindex];
                if (left_edge == -1)
                    {
                        //Variant 0 can only break homozygosity
                        //run if it is an xton
                        if (xtons.front() == 0 && hapi[0] != hapj[0])
                            {
                                left_edge = 0;
                                edges[lindex] = left_edge;
                                rv = true;
                            }
                    }
                if (left_edge >= 0)
                    {
                        std::size_t rindex = i * m.nsam + j;
                        std::int64_t right_edge = edges[rindex];
                        if (right_edge == -1
                            || static_cast<std::size_t>(right_edge) <= core)
                            {
                                rv = false;
                                // To update right edge:
                                // Iterate over all xtons > core
                                // and check if haplotypes i,j
                                // differ at those sites.
                                for (auto r = flanks.second;
                                     r < xtons.size() && rv == false; ++r)
                                    {
                                        if (!(core <= xtons[r]))
                                            {
                                                throw std::logic_error(
                                                    "r < core");
                                            }
                                        if (hapi[xtons[r]] != hapj[xtons[r]]
                                            && !(hapi[xtons[r]] < 0))
                                            {
                                                rv = true;
                                                edges[rindex] = r;
                                            }
                                    }
                            }
                    }
            }
        else if (std::binary_search(xtons.begin() + flanks.first,
                                    xtons.begin() + flanks.second, core))
            {
                //seqs i and j differ and core is an xton,
                //thus core is a new left
                edges[lindex] = core;
            }
        return rv;
    }
    std::pair<std::int64_t, std::int64_t>
    get_flanking_xtons(const std::vector<int64_t>& xtons,
                       const std::size_t core)
    {
        auto left_xton = std::upper_bound(
            xtons.rbegin(), xtons.rend(), core,
            [](const std::int64_t element, const std::size_t value) {
                return element > value;
            });
        if (left_xton == xtons.rend())
            {
                return std::make_pair(-1, -1);
            }
        auto right_xton = std::upper_bound(
            left_xton.base(), xtons.end(), core,
            [](const std::int64_t element, const std::size_t value) {
                return element < value;
            });
        if (right_xton == xtons.end())
            {
                return std::make_pair(-1, -1);
            }
        return std::make_pair(std::distance(xtons.begin(), left_xton.base()),
                              std::distance(xtons.begin(), right_xton));
    }
} // namespace

namespace Sequence
{
    //TODO refactor to use suffix_edges
    std::vector<nSLiHS>
    nslx(const VariantMatrix& m, const std::int8_t refstate, const int x)
    {
        //Need to get indexes of all x-tons.
        //Then, if two seqs differ at an x-ton,
        //the stats get updated.
        std::vector<std::int64_t> xtons;
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(m.nsites); ++i)
            {
                auto r = get_ConstRowView(m, static_cast<std::size_t>(i));
                auto nonref = std::count_if(
                    r.begin(), r.end(), [refstate](const std::int8_t a) {
                        return a != refstate && !(a < 0);
                    });
                if (nonref == x)
                    {
                        xtons.push_back(i);
                    }
            }
        std::vector<nSLiHS> rv;
        if (xtons.empty() || !m.nsam || !m.nsites)
            {
                return rv;
            }

        std::vector<std::int64_t> edges(m.nsam * m.nsam, -1);
        std::size_t lindex, rindex;
        for (std::size_t core = 0; core < m.nsites; ++core)
            {
                auto core_view = get_ConstRowView(m, core);
                // Doing any work requires the existence
                // of x-tons left and right of core
                auto flanks = get_flanking_xtons(xtons, core);
                double nsl_values[2] = { 0, 0 };
                double ihs_values[2] = { 0, 0 };
                int counts[2] = { 0, 0 };
                if (flanks.first != -1)
                    {
                        for (std::size_t i = 0; i < m.nsam - 1; ++i)
                            {
                                auto sample_i = get_ConstColView(m, i);
                                for (std::size_t j = i + 1; j < m.nsam; ++j)
                                    {
                                        if (update_edge_matrix(
                                                m, xtons, flanks, core_view,
                                                sample_i, edges, core, i, j))
                                            {
                                                lindex = j * m.nsam + i;
                                                rindex = i * m.nsam + j;
                                                summstats_details::update_counts(
                                                    nsl_values, ihs_values,
                                                    counts, m.nsites,
                                                    m.positions,
                                                    static_cast<std::size_t>(
                                                        core_view[i]
                                                        == refstate),
                                                    edges[lindex],
                                                    edges[rindex]);
                                            }
                                    }
                            }
                    }
                rv.emplace_back(summstats_details::get_stat(core_view, refstate, nsl_values,
                                         ihs_values, counts));
            }
        return rv;
    }
} // namespace Sequence
