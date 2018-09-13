//! \file testGarudStatistics.cc @brief unit tests for H1, H12, H2/H1 stats

#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/garud.hpp>
#include <Sequence/summstats/classics.hpp>
#include "VariantMatrixFixture.hpp"
#include <boost/test/unit_test.hpp>

BOOST_FIXTURE_TEST_SUITE(test_garud_stats, dataset)

BOOST_AUTO_TEST_CASE(test_garud_stats)
/* This test compares the results to alternative
 * (and less efficient) implementations
 * of the same procedures.
 */
{
    auto G = Sequence::garud_statistics(m);

    //Their H1 stat is the same as 1 - haplotype diversity.
    auto hdiv = Sequence::haplotype_diversity(m);
    // This test is a cheat, as this is exactly
    // how H1 is calculated internally
    BOOST_REQUIRE_EQUAL(1.0 - hdiv, G.H1);

    //Let's do an alternative method.
    std::vector<std::string> haps;
    for (auto i = 0; i < m.nsam; ++i)
        {
            auto c = Sequence::get_ConstColView(m, i);
            std::string h;
            for (auto ci : c)
                {
                    h += (ci == 0) ? '0' : '1';
                }
            haps.push_back(std::move(h));
        }
    // unique haplotypes
    std::set<std::string> uhaps(haps.begin(), haps.end());
    std::vector<int> hcounts;
    for (auto& u : uhaps)
        {
            hcounts.push_back(std::count(haps.begin(), haps.end(), u));
        }
    std::sort(hcounts.begin(), hcounts.end(), std::greater<int>());
    double H1 = 0.0;
    double nsam = static_cast<double>(m.nsam);
    for (auto hc : hcounts)
        {
            H1 += static_cast<double>(hc) * static_cast<double>(hc - 1);
        }
    H1 /= ((nsam) * (nsam - 1));

    //We cannot require equality because we are summing things in
    //different orders.  Instead, we check that relative error is
    //very low.
    BOOST_CHECK_CLOSE(H1, G.H1, 1e-6);

    //H12 is H1, but combining the freq of the two most common haplotypes
    decltype(hcounts) hcounts2(hcounts.begin() + 1, hcounts.end());
    hcounts2[0] += hcounts[0];
    double H12 = 0.0;
    for (auto hc : hcounts2)
        {
            double x = static_cast<double>(hc) * static_cast<double>(hc - 1);
            x /= (nsam * (nsam - 1));
            H12 += x;
        }
    BOOST_CHECK_CLOSE(H12, G.H12, 1e-6);

    double p1sq = static_cast<double>(hcounts[0] * (hcounts[0] - 1))
                  / (nsam * (nsam - 1.0));
    double H2 = H1 - p1sq;
    BOOST_CHECK_CLOSE(H2 / H1, G.H2H1, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
