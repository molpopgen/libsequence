//! \file testLD.cc @brief unit tests for LD-related calculations

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/ld.hpp>
#include "VariantMatrixFixture.hpp"
#include <boost/test/unit_test.hpp>

BOOST_FIXTURE_TEST_SUITE(test_LD, dataset)

BOOST_AUTO_TEST_CASE(test_two_locus_haplotype_counts)
{
    std::vector<int> results;
    for (std::size_t i = 0; i < m.nsites - 1; ++i)
        {
            for (std::size_t j = i + 1; j < m.nsites; ++j)
                {
                    std::vector<std::pair<std::int8_t, std::int8_t>> haps;
                    auto hc
                        = Sequence::two_locus_haplotype_counts(m, i, j, true);
                    auto ri = Sequence::get_ConstRowView(m, i);
                    auto rj = Sequence::get_ConstRowView(m, j);
                    for (std::size_t k = 0; k < ri.size(); ++k)
                        {
                            haps.emplace_back(ri[k], rj[k]);
                        }
                    std::sort(haps.begin(), haps.end());
                    auto end_of_unique_haps
                        = std::unique(haps.begin(), haps.end());
                    BOOST_REQUIRE_EQUAL(
                        hc.size(), static_cast<std::size_t>(std::distance(
                                       haps.begin(), end_of_unique_haps)));
                }
        }
}

BOOST_AUTO_TEST_SUITE_END()
