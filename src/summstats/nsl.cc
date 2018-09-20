#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/nsl.hpp>
#include "nsl_common.hpp"
#include "algorithm.hpp"

/// \example nSL_from_ms.cc
namespace
{
    std::int64_t
    get_left(const Sequence::ConstColView& sample_i,
             const Sequence::ConstColView& sample_j, const std::size_t core,
             std::int64_t minleft)
    {
        std::int64_t left = static_cast<std::int64_t>(core) - 1;
        auto state_i = sample_i.begin() + left;
        auto state_j = sample_j.begin() + left;
        while (left >= minleft)
            {
                if (*state_i != *state_j && !(*state_i < 0 || *state_j < 0))
                    {
                        break;
                    }
                --left;
                --state_i;
                --state_j;
            }
        return left;
    }

    std::int64_t
    get_right(const Sequence::ConstColView& sample_i,
              const Sequence::ConstColView& sample_j, const std::size_t core)
    {
        auto state_i = sample_i.begin() + static_cast<std::int64_t>(core) + 1;
        auto state_j = sample_j.begin() + static_cast<std::int64_t>(core) + 1;
        auto m = Sequence::summstats_algo::mismatch_skip_missing(
            state_i, sample_i.end(), state_j);
        return std::distance(sample_i.begin(), m.first);
    }

    inline bool
    update_edge_matrix(const Sequence::VariantMatrix& m,
                       const Sequence::ConstRowView& core_view,
                       const Sequence::ConstColView& hapi,
                       const Sequence::ConstColView& hapj,
                       Sequence::summstats_details::suffix_edges& edges,
                       const std::size_t core, const std::size_t i,
                       const std::size_t j)
    // For the nsl operations, this function keeps track
    // of the left/right boundaries of haplotype homozygosity
    // as core moves through the sample:
    {
        bool rv = false;
        if (core_view[i] == core_view[j] && !(core_view[i] < 0))
            {
                rv = true;
                if (edges.left >= 0)
                    {
                        // This next line is a gotcha
                        // to to signed/unsigned comparison:
                        if (edges.right == -1 || edges.right <= core)
                            {
                                edges.right = get_right(hapi, hapj, core);
                            }
                    }
            }
        else if (!(core_view[i] < 0 || core_view[j] < 0))
            {
                edges.left = static_cast<std::int64_t>(core);
            }
        return rv;
    }
} // namespace

namespace Sequence
{
    nSLiHS
    nsl(const VariantMatrix& m, const std::size_t core,
        const std::int8_t refstate)
    {
        auto core_view = get_ConstRowView(m, core);
        // Keep track of distances from core site
        // for nsl and ihs separately
        double nsl_values[2] = { 0, 0 };
        double ihs_values[2] = { 0, 0 };
        //Count sample size for non-ref and ref alleles at core site contributing to nSL
        int counts[2] = { 0, 0 };
        for (std::size_t i = 0; i < m.nsam - 1; ++i)
            {
                auto sample_i = get_ConstColView(m, i);
                for (std::size_t j = i + 1; j < m.nsam; ++j)
                    {
                        if (core_view[i] == core_view[j] && core_view[i] >= 0)
                            {
                                auto sample_j = get_ConstColView(m, j);
                                //Find where samples i and j differ
                                auto left
                                    = get_left(sample_i, sample_j, core, 0);
                                if (left >= 0)
                                    {
                                        auto right = get_right(sample_i,
                                                               sample_j, core);
                                        summstats_details::update_counts(
                                            nsl_values, ihs_values, counts,
                                            m.nsites, m.positions,
                                            static_cast<std::size_t>(
                                                core_view[i] == refstate),
                                            left, right);
                                    }
                            }
                    }
            }
        return summstats_details::get_stat(core_view, refstate, nsl_values,
                                           ihs_values, counts);
    }

    std::vector<nSLiHS>
    nsl(const VariantMatrix& m, const std::int8_t refstate)
    {
        std::vector<nSLiHS> rv;
        if (m.nsam == 0)
            {
                return rv;
            }
        rv.reserve(m.nsites);

        // A matrix keeping track of the
        // index where sample i,j last differed.
        // The lower left corresponds to left edges,
        // and the upper right are the right edges.
        // -1 mean unevaluated.
        auto npairs = m.nsam * (m.nsam - 1) / 2;
        std::vector<summstats_details::suffix_edges> edges(npairs);
        std::size_t lindex, rindex;
        std::vector<ConstColView> alleles;
        alleles.reserve(m.nsam);
        for (std::size_t i = 0; i < m.nsam; ++i)
            {
                alleles.push_back(get_ConstColView(m, i));
            }
        for (std::size_t core = 0; core < m.nsites; ++core)
            {
                std::size_t pair_index = 0;
                auto core_view = get_ConstRowView(m, core);
                double nsl_values[2] = { 0, 0 };
                double ihs_values[2] = { 0, 0 };
                //Count sample size for non-ref and
                //ref alleles at core site contributing to nSL
                int counts[2] = { 0, 0 };
                for (std::size_t i = 0; i < m.nsam - 1; ++i)
                    {
                        const auto& hapi = alleles[i];
                        for (std::size_t j = i + 1; j < m.nsam;
                             ++j, ++pair_index)
                            {
                                if (update_edge_matrix(
                                        m, core_view, hapi, alleles[j],
                                        edges[pair_index], core, i, j))
                                    {
                                        lindex = j * m.nsam + i;
                                        rindex = i * m.nsam + j;
                                        summstats_details::update_counts(
                                            nsl_values, ihs_values, counts,
                                            m.nsites, m.positions,
                                            static_cast<std::size_t>(
                                                core_view[i] == refstate),
                                            edges[pair_index].left,
                                            edges[pair_index].right);
                                    }
                            }
                    }
                rv.emplace_back(summstats_details::get_stat(
                    core_view, refstate, nsl_values, ihs_values, counts));
            }
        return rv;
    }
} // namespace Sequence
