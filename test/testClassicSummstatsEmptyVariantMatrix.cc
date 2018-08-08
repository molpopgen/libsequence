/*! \file testClassicSummstatsEmptyVariantMatrix.cc @brief unit tests for pop gen summary statistics
 *
 *  When a VariantMatrix is empty, there are two possibilities:
 *
 *  1. There are no variable sites.  Equivalently, all variant sites
 *     may have been filtered out for some reason, resulting in
 *     an empty matrix
 *  2. There are sites in the matrix, but they are all invariant (modulo missing data).
 *  2. The sample size is zero.
 *
 *  It is not possible to distinguish the three cases from looking
 *  at the contents of a VariantMatrix.  Thus, libsequence assumes
 *  that case 1 or 2 applies.
 *
 */

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/classics.hpp>
#include "VariantMatrixFixture.hpp"
#include <boost/test/unit_test.hpp>

BOOST_FIXTURE_TEST_SUITE(test_classic_stats_with_empty_variant_matrix,
                         invariantdataset)

BOOST_AUTO_TEST_CASE(test_thetapi)
{
    auto pi = Sequence::thetapi(empty_counts);
    BOOST_REQUIRE_EQUAL(pi, 0.0);
    pi = Sequence::thetapi(invariant_counts);
    BOOST_REQUIRE_EQUAL(pi, 0.0);
}

BOOST_AUTO_TEST_CASE(test_thetaw)
{
    auto tw = Sequence::thetaw(empty_counts);
    BOOST_REQUIRE_EQUAL(tw, 0.0);
    tw = Sequence::thetaw(invariant_counts);
    BOOST_REQUIRE_EQUAL(tw, 0.0);
}

BOOST_AUTO_TEST_CASE(test_thetah)
{
    auto th = Sequence::thetah(empty_counts, 0);
    BOOST_REQUIRE_EQUAL(th, 0.0);
    th = Sequence::thetah(invariant_counts, 0);
    BOOST_REQUIRE_EQUAL(th, 0.0);
}

BOOST_AUTO_TEST_CASE(test_thetal)
{
    auto tl = Sequence::thetal(empty_counts, 0);
    BOOST_REQUIRE_EQUAL(tl, 0.0);
    tl = Sequence::thetal(invariant_counts, 0);
    BOOST_REQUIRE_EQUAL(tl, 0.0);
}

BOOST_AUTO_TEST_CASE(test_tajd)
{
    auto td = Sequence::tajd(empty_counts);
    BOOST_REQUIRE_EQUAL(std::isnan(td), true);
    td = Sequence::tajd(invariant_counts);
    BOOST_REQUIRE_EQUAL(std::isnan(td), true);
}

BOOST_AUTO_TEST_CASE(test_faywuh)
{
    auto fwh = Sequence::faywuh(empty_counts, 0);
    BOOST_REQUIRE_EQUAL(std::isnan(fwh), true);
    fwh = Sequence::faywuh(invariant_counts, 0);
    BOOST_REQUIRE_EQUAL(std::isnan(fwh), true);
}

BOOST_AUTO_TEST_CASE(test_hprime)
{
    auto hp = Sequence::hprime(empty_counts, 0);
    BOOST_REQUIRE_EQUAL(std::isnan(hp), true);
    hp = Sequence::hprime(invariant_counts, 0);
    BOOST_REQUIRE_EQUAL(std::isnan(hp), true);
}

// Tests of haplotype statistics are more complex.
// It is not as obvious (to me) what to do for an empty
// matrix.

BOOST_AUTO_TEST_CASE(test_num_haplotypes)
{
    auto nh = Sequence::number_of_haplotypes(empty);
    BOOST_REQUIRE_EQUAL(nh, -1);
    nh = Sequence::number_of_haplotypes(invariant);
    BOOST_REQUIRE_EQUAL(nh, 1);
}

BOOST_AUTO_TEST_CASE(test_haplotype_diversity)
{
    auto hd = Sequence::haplotype_diversity(empty);
    BOOST_REQUIRE_EQUAL(std::isnan(hd), true);
    hd = Sequence::haplotype_diversity(invariant);
    BOOST_REQUIRE_EQUAL(hd, 0.0);
}

BOOST_AUTO_TEST_CASE(test_rmin)
{
    auto rm = Sequence::rmin(empty);
    BOOST_CHECK_EQUAL(rm, -1);
    rm = Sequence::rmin(invariant);
    BOOST_CHECK_EQUAL(rm, 0);
}

BOOST_AUTO_TEST_SUITE_END()
