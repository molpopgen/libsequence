#include <algorithm>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <Sequence/summstats/util.hpp>
#include <Sequence/summstats/generic.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include "algorithm.hpp"

namespace Sequence
{
    //TODO: modify to ignore sequences
    //with missing data above a certain threshold?
	//TODO: needs unit test
    std::vector<std::int32_t>
    difference_matrix(const VariantMatrix& m)
    {
        std::vector<std::int32_t> rv;
        rv.reserve(m.nsam);
        for (std::size_t i = 0; i < m.nsam - 1; ++i)
            {
                auto ci = get_ConstColView(m, i);
                for (std::size_t j = i + 1; j < m.nsam; ++j)
                    {
                        auto cj = get_ConstColView(m, j);
                        std::int32_t ndiffs
                            = summstats_algo::ndiff_skip_missing(
                                ci.begin(), ci.end(), cj.begin());
                        rv.push_back(ndiffs);
                    }
            }
        return rv;
    }

    std::vector<std::int32_t>
    is_different_matrix(const VariantMatrix& m)
    {
        std::vector<std::int32_t> rv;
        rv.reserve(m.nsam);
        std::vector<ConstColView> alleles;
        alleles.reserve(m.nsam);
        for (std::size_t i = 0; i < m.nsam; ++i)
            {
                alleles.push_back(get_ConstColView(m, i));
            }
        for (std::size_t i = 0; i < m.nsam - 1; ++i)
            {
                for (std::size_t j = i + 1; j < m.nsam; ++j)
                    {
                        auto m = summstats_algo::mismatch_skip_missing(
                            alleles[i].begin(), alleles[i].end(),
                            alleles[j].begin());
                        rv.push_back(m.first != alleles[i].end());
                    }
            }
        return rv;
    }

    std::vector<std::int32_t>
    label_haplotypes(const VariantMatrix& m)
    {
        std::vector<std::int32_t> rv(m.nsam, -1);
        if (rv.empty())
            {
                return rv;
            }
        rv.reserve(m.nsam);
        const auto dm = is_different_matrix(m);
        auto dmi = dm.cbegin();
        int next_label = 0;
        // We got all the way to nsam for the
        // case where the last haplotype is unique.
        // In 1.9.4, we went to nsam-1, which was wrong.
        // Fixed in 1.9.5
        for (std::size_t i = 0; i < m.nsam; ++i)
            {
                if (rv[i] < 0)
                    {
                        if (!all_missing(get_ConstColView(m, i)))
                            {
                                rv[i] = next_label;
                                for (std::size_t j = i + 1; j < m.nsam;
                                     ++j, ++dmi)
                                    {
                                        if (!all_missing(
                                                get_ConstColView(m, j)))
                                            {
                                                if (*dmi == 0)
                                                    {
                                                        rv[j] = next_label;
                                                    }
                                            }
                                    }
                                ++next_label;
                            }
                    }
            }
        return rv;
    }

    std::int32_t
    number_of_haplotypes(const VariantMatrix& m)
    {
        if (m.data.empty() || !m.nsam)
            {
                return -1;
            }
        auto labels = label_haplotypes(m);
        labels.erase(
            std::remove_if(labels.begin(), labels.end(),
                           [](decltype(labels[0]) x) { return x < 0; }),
            labels.end());
        std::sort(labels.begin(), labels.end());
        auto u = std::unique(labels.begin(), labels.end());
        return static_cast<std::int32_t>(std::distance(labels.begin(), u));
    }

    double
    haplotype_diversity(const VariantMatrix& m)
    {
        if (m.data.empty() || !m.nsam)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto labels = label_haplotypes(m);
        auto nmissing
            = std::count_if(labels.begin(), labels.end(),
                            [](decltype(labels[0]) x) { return x < 0; });
        auto nsam_adjusted = m.nsam - static_cast<decltype(m.nsam)>(nmissing);
        labels.erase(
            std::remove_if(labels.begin(), labels.end(),
                           [](decltype(labels[0]) x) { return x < 0; }),
            labels.end());
        std::unordered_map<std::int32_t, std::int32_t> label_counts;
        for (auto l : labels)
            {
                label_counts[l]++;
            }
        return diversity_from_counts(label_counts, nsam_adjusted);
    }

} // namespace Sequence
