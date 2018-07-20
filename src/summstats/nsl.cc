#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/nsl.hpp>

namespace
{
    std::int64_t
    get_left(const Sequence::ConstColView& sample_i,
             const Sequence::ConstColView& sample_j, const std::size_t core)
    {
        std::int64_t left = static_cast<std::int64_t>(core) - 1;
        std::size_t left_index;

        while (left >= 0)
            {
                left_index = static_cast<std::size_t>(left);
                if (sample_i[left_index] >= 0 && sample_j[left_index] >= 0
                    && sample_i[left_index] != sample_j[left_index])
                    {
                        break;
                    }
                --left;
            }
        return left;
    }

    std::int64_t
    get_right(const Sequence::ConstColView& sample_i,
              const Sequence::ConstColView& sample_j, const std::size_t core,
              const std::int64_t nsites)
    {
        std::int64_t right = static_cast<std::int64_t>(core) + 1;
        std::size_t right_index;
        while (right < nsites)
            {
                right_index = static_cast<std::size_t>(right);
                if (sample_i[right_index] >= 0 && sample_j[right_index] >= 0
                    && sample_i[right_index] != sample_j[right_index])
                    {
                        break;
                    }
                ++right;
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
                    += positions[static_cast<std::size_t>(left)]
                       - positions[static_cast<std::size_t>(right)];
                counts[index]++;
            }
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
                        if (core_view[i] == core_view[j])
                            {
                                auto sample_j = get_ConstColView(m, j);
                                //Find where samples i and j differ
                                auto left = get_left(sample_i, sample_j, core);
                                auto right = get_right(
                                    sample_i, sample_j, core,
                                    static_cast<std::int64_t>(m.nsites));
                                update_counts(nsl_values, ihs_values, counts,
                                              m.nsites, m.positions,
                                              static_cast<std::size_t>(
                                                  core_view[i] == refstate),
                                              left, right);
                            }
                    }
            }
        double nSL_den = nsl_values[0] / static_cast<double>(counts[0]);
        double nSL_num = nsl_values[1] / static_cast<double>(counts[1]);
        double iHS_den = nsl_values[0] / static_cast<double>(counts[0]);
        double iHS_num = nsl_values[1] / static_cast<double>(counts[1]);
        auto nonrefcount = static_cast<std::int32_t>(
            std::count_if(core_view.begin(), core_view.end(),
                          [refstate](const std::int8_t i) {
                              return i >= 0 && i != refstate;
                          }));
        return nSLiHS{ std::log(nSL_num) - std::log(nSL_den),
                       std::log(iHS_num) - std::log(iHS_den), nonrefcount };
    }
} // namespace Sequence
