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
                    auto hc
                        = Sequence::two_locus_haplotype_counts(m, i, j, true);
                    results.push_back(hc.size());
                }
        }
    BOOST_REQUIRE_EQUAL(results[0], 4);
    BOOST_REQUIRE_EQUAL(results[1], 4);
    BOOST_REQUIRE_EQUAL(results[2], 3);
}

BOOST_AUTO_TEST_SUITE_END()
