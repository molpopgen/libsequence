#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/nsl.hpp>

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
              const Sequence::ConstColView& sample_j, const std::size_t core,
              const std::int64_t nsites)
    {
        std::int64_t right = static_cast<std::int64_t>(core) + 1;
        auto state_i = sample_i.begin() + right;
        auto state_j = sample_j.begin() + right;
        while (right < nsites)
            {
                if (*state_i != *state_j && !(*state_i < 0 || *state_j < 0))
                    {
                        break;
                    }
                ++right;
                ++state_i;
                ++state_j;
            }
        return right;
    }

    void
    update_counts(double nsl_values[2], double ihs_values[2], int counts[2],
                  const std::size_t nsites,
                  const std::vector<double>& positions,
                  const std::size_t index, const std::int64_t left,
                  const std::int64_t right)
    {
        if (left >= 0 && static_cast<std::size_t>(right) < nsites)
            //Then there are SNPs differentiating
            //i and j within the region
            {
                nsl_values[index] += static_cast<double>(right - left);
                //TODO: check if we need to add one?
                ihs_values[index]
                    += positions[static_cast<std::size_t>(right)]
                       - positions[static_cast<std::size_t>(left)];
                counts[index]++;
            }
    }

    inline void
    update_edge_matrix(const Sequence::VariantMatrix& m,
                       const Sequence::ConstColView& hapi,
                       const Sequence::ConstColView& hapj,
                       std::vector<std::int64_t>& edges,
                       const std::size_t core, const std::size_t i,
                       const std::size_t j)
    // For the all_nsl operations, this function keeps track
    // of the left/right boundaries of haplotype homozygosity
    // as core moves through the sample:
    {
        std::size_t lindex = j * m.nsam + i;
        std::int64_t left_edge = edges[lindex];
        if (left_edge == -1)
            {
                left_edge = get_left(hapi, hapj, core, 0);
                edges[lindex] = left_edge;
            }
        if (left_edge >= 0)
            {
                std::size_t rindex = i * m.nsam + j;
                std::int64_t right_edge = edges[rindex];
                // This next line is a gotcha
                // to to signed/unsigned comparison:
                if (right_edge == -1
                    || static_cast<std::size_t>(right_edge) <= core)
                    {
                        right_edge
                            = get_right(hapi, hapj, core,
                                        static_cast<std::int64_t>(m.nsites));
                        edges[rindex] = right_edge;
                    }
            }
    }

    inline Sequence::nSLiHS
    get_stat(const Sequence::ConstRowView& core_view,
             const std::int8_t refstate, const double nsl_values[2],
             const double ihs_values[2], const int counts[2])
    {

        double nSL_den = nsl_values[0] / static_cast<double>(counts[0]);
        double nSL_num = nsl_values[1] / static_cast<double>(counts[1]);
        double iHS_den = ihs_values[0] / static_cast<double>(counts[0]);
        double iHS_num = ihs_values[1] / static_cast<double>(counts[1]);
        auto nonrefcount = static_cast<std::int32_t>(
            std::count_if(core_view.begin(), core_view.end(),
                          [refstate](const std::int8_t i) {
                              return i >= 0 && i != refstate;
                          }));
        return Sequence::nSLiHS{ std::log(nSL_num) - std::log(nSL_den),
                                 std::log(iHS_num) - std::log(iHS_den),
                                 nonrefcount };
    }

} // namespace

namespace Sequence
{
    nSLiHS
    nsl(const VariantMatrix& m, const std::size_t core,
        const std::int8_t refstate)
    /*! \brief nSL and iHS statistics
     * \param m A VariantMatrix
     * \param core The index of the core site
     * \param refstate The value of the reference/ancestral allelic state
     *
     * \return an nSLiHS object
     * \ingroup popgenanalysis
     *
     * See nSL_from_ms.cc for example
     *
     * See \cite Ferrer-Admetlla2014-wa for details.
     */
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
                                        auto right = get_right(
                                            sample_i, sample_j, core,
                                            static_cast<std::int64_t>(
                                                m.nsites));
                                        update_counts(
                                            nsl_values, ihs_values, counts,
                                            m.nsites, m.positions,
                                            static_cast<std::size_t>(
                                                core_view[i] == refstate),
                                            left, right);
                                    }
                            }
                    }
            }
        return get_stat(core_view, refstate, nsl_values, ihs_values, counts);
    }

    // TODO: A further optimization may be possible.
    // If you move the iteration over cores to the innermost
    // loop, then you can keep recording the same value
    // up until you hit the right position.  The
    // reason is that you know that homozygosity
    // extends to the right and all you have to
    // check for is ancestral, derived, or missing.
    // This change would reduce the extra n^2 memory
    // required but it would introduce an O(nsites)
    // traversal for all sequence pairs.
    std::vector<nSLiHS>
    all_nsl(const VariantMatrix& m, const std::int8_t refstate)
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
        std::vector<std::int64_t> edges(m.nsam * m.nsam, -1);
        std::size_t lindex, rindex;
        for (std::size_t core = 0; core < m.nsites; ++core)
            {
                auto core_view = get_ConstRowView(m, core);
                double nsl_values[2] = { 0, 0 };
                double ihs_values[2] = { 0, 0 };
                //Count sample size for non-ref and
                //ref alleles at core site contributing to nSL
                int counts[2] = { 0, 0 };
                for (std::size_t i = 0; i < m.nsam - 1; ++i)
                    {
                        auto sample_i = get_ConstColView(m, i);
                        for (std::size_t j = i + 1; j < m.nsam; ++j)
                            {
                                if (core_view[i] == core_view[j]
                                    && core_view[i] >= 0)
                                    {
                                        auto sample_j = get_ConstColView(m, j);
                                        lindex = j * m.nsam + i;
                                        rindex = i * m.nsam + j;
                                        update_edge_matrix(m, sample_i,
                                                           sample_j, edges,
                                                           core, i, j);
                                        update_counts(
                                            nsl_values, ihs_values, counts,
                                            m.nsites, m.positions,
                                            static_cast<std::size_t>(
                                                core_view[i] == refstate),
                                            edges[lindex], edges[rindex]);
                                    }
                                else //The sites differ, making this a new left edge
                                    {
                                        edges[j * m.nsam + i]
                                            = static_cast<std::int64_t>(core);
                                    }
                            }
                    }
                rv.emplace_back(get_stat(core_view, refstate, nsl_values,
                                         ihs_values, counts));
            }
        return rv;
    }

    nSLiHS
    nslx(const VariantMatrix& m, const std::size_t core,
         const std::int8_t refstate, const int x)
    {
        //Need to get indexes of all x-tons.
        //Then, if two seqs differ at an x-ton,
        //the stats get updated.
    }
} // namespace Sequence
