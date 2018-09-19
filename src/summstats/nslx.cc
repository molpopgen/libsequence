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
                       const Sequence::ConstColView& hapj,
                       Sequence::summstats_details::suffix_edges& edges,
                       const std::size_t core, const std::size_t i,
                       const std::size_t j)
    // edge updating for nslx
    // precondition: !xtons.empty()
    {
        bool rv = false;
        if (core_view[i] == core_view[j] && !(core_view[i] < 0))
            {
                if (edges.left == -1)
                    {
                        //Variant 0 can only break homozygosity
                        //run if it is an xton
                        if (xtons.front() == 0 && hapi[0] != hapj[0]
                            && !(hapi[0] < 0 || hapj[0] < 0))
                            {
                                edges.left = 0;
                                rv = true;
                            }
                    }
                if (edges.left >= 0)
                    {
                        if (edges.right == -1
                            || static_cast<std::size_t>(edges.right) <= core)
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
                                            && !(hapi[xtons[r]] < 0
                                                 || hapj[xtons[r]] < 0))
                                            {
                                                rv = true;
                                                edges.right = r;
                                            }
                                    }
                            }
                    }
            }
        else if (!(hapi[core] < 0 || hapj[core] < 0)
                 && std::binary_search(xtons.begin() + flanks.first,
                                       xtons.begin() + flanks.second, core))
            {
                //seqs i and j differ and core is an xton,
                //thus core is a new left
                edges.left = core;
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

        std::size_t npairs = m.nsam * (m.nsam - 1) / 2;
        std::vector<summstats_details::suffix_edges> edges(npairs);
        std::vector<ConstColView> alleles;
        alleles.reserve(m.nsam);
        for (std::size_t i = 0; i < m.nsam; ++i)
            {
                alleles.push_back(get_ColView(m, i));
            }
        for (std::size_t core = 0; core < m.nsites; ++core)
            {
                auto core_view = get_ConstRowView(m, core);
                // Doing any work requires the existence
                // of x-tons left and right of core
                auto flanks = get_flanking_xtons(xtons, core);
                double nsl_values[2] = { 0, 0 };
                double ihs_values[2] = { 0, 0 };
                int counts[2] = { 0, 0 };
                std::size_t pair_index = 0;
                if (flanks.first != -1)
                    {
                        for (std::size_t i = 0; i < m.nsam - 1; ++i)
                            {
                                for (std::size_t j = i + 1; j < m.nsam;
                                     ++j, ++pair_index)
                                    {
                                        if (update_edge_matrix(
                                                m, xtons, flanks, core_view,
                                                alleles[i], alleles[j],
                                                edges[pair_index],
                                                core, i, j))
                                            {
                                                summstats_details::
                                                    update_counts(
                                                        nsl_values, ihs_values,
                                                        counts, m.nsites,
                                                        m.positions,
                                                        static_cast<
                                                            std::size_t>(
                                                            core_view[i]
                                                            == refstate),
                                                        edges[pair_index].left,
                                                        edges[pair_index]
                                                            .right);
                                            }
                                    }
                            }
                    }
                rv.emplace_back(summstats_details::get_stat(
                    core_view, refstate, nsl_values, ihs_values, counts));
            }
        return rv;
    }
} // namespace Sequence
