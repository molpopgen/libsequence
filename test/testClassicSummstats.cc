//! \file testClassicSummstats.cc @brief unit tests for pop gen summary statistics

#include <cmath>
#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/classics.hpp>
#include "VariantMatrixFixture.hpp"
#include <boost/test/unit_test.hpp>

static double
manual_pi(const Sequence::VariantMatrix& m)
{
    double manual = 0.0;
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto r = Sequence::get_ConstRowView(m, i);
            int ndiffs = 0, ncomps = 0;
            for (std::size_t j = 0; j < m.nsam - 1; ++j)
                {
                    for (std::size_t k = j + 1; k < m.nsam; ++k)
                        {
                            if (r[j] >= 0 && r[k] >= 0)
                                {
                                    if (r[j] != r[k])
                                        {
                                            ++ndiffs;
                                        }
                                    ++ncomps;
                                }
                        }
                }
            manual
                += static_cast<double>(ndiffs) / static_cast<double>(ncomps);
        }
    return manual;
}

static double
manual_thetah(const Sequence::VariantMatrix& m, const std::int8_t refstate)
{
    double nnm1 = static_cast<double>(m.nsam * (m.nsam - 1));
    double h = 0.0;
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto r = Sequence::get_ConstRowView(m, i);
            unsigned nnonref = 0;
            for (auto ri : r)
                {
                    if (ri != refstate)
                        {
                            ++nnonref;
                        }
                }
            h += std::pow(nnonref, 2.0);
        }
    return h / nnm1;
}

BOOST_FIXTURE_TEST_SUITE(test_classic_stats, dataset)

BOOST_AUTO_TEST_CASE(test_thetapi)
{
    auto pi = Sequence::thetapi(c);
    auto manual = manual_pi(m);
    // Cannot require equal b/c we aren't doing ops
    // in same order.
    BOOST_CHECK_CLOSE(pi, manual, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_thetapi_with_mising_data)
// If this test passes, it implies
// that Sequence::StateCounter behaves
// correctly w.r.to handling missing data.
{
    // Make a bunch of missing data, all w/
    // different missing data encoding
    for (int i = 1; i < static_cast<int>(m.nsites); i += 2)
        {
            m.data[i] = -i;
        }
    auto pi = Sequence::thetapi(Sequence::AlleleCountMatrix(m));

    auto manual = manual_pi(m);
    // Cannot require equal b/c we aren't doing ops
    // in same order.
    BOOST_CHECK_CLOSE(pi, manual, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_num_poly_sites)
{
    auto S = Sequence::nvariable_sites(m);
    // Not universally true, but is true here:
    BOOST_REQUIRE_EQUAL(m.nsites, S);

    auto tm = Sequence::total_number_of_mutations(m);
    // Not universally true, but is true here:
    BOOST_REQUIRE_EQUAL(m.nsites, tm);

    //Make the last site invariant:
    auto r = Sequence::get_RowView(m, m.nsites - 1);
    std::fill(r.begin(), r.end(), 0);
    S = Sequence::nvariable_sites(m);
    BOOST_REQUIRE_EQUAL(S, m.nsites - 1);
    tm = Sequence::total_number_of_mutations(m);
    BOOST_REQUIRE_EQUAL(m.nsites - 1, tm);
}

BOOST_AUTO_TEST_CASE(test_total_num_mutations)
{
    auto r = Sequence::get_RowView(m, m.nsites - 1);
    // Add a third character state
    r[2] = 2;
    auto tm = Sequence::total_number_of_mutations(m);
    BOOST_REQUIRE_EQUAL(tm, m.nsites + 1);
}

BOOST_AUTO_TEST_CASE(test_nbiallelic_sites)
{
    auto S2 = Sequence::nbiallelic_sites(m);
    BOOST_REQUIRE_EQUAL(S2, m.nsites);
    auto r = Sequence::get_RowView(m, m.nsites - 1);
    // Add a third character state
    r[2] = 2;
    S2 = Sequence::nbiallelic_sites(m);
    BOOST_REQUIRE_EQUAL(S2, m.nsites - 1);
}

BOOST_AUTO_TEST_CASE(test_count_alleles)
{
    auto S2 = Sequence::nbiallelic_sites(m);
    auto ac = Sequence::allele_counts(m);
    auto nb = 0;
    for (auto i : ac)
        {
            if (i.nstates == 2)
                {
                    ++nb;
                }
        }
    BOOST_REQUIRE_EQUAL(S2, nb);
}

BOOST_AUTO_TEST_CASE(test_thetaw)
// The above two tests imply that
// thetaw is correct.  Thus,
// all we're doing below is
// checking the denominator.
{
    auto w = Sequence::thetaw(c);
    double S = m.nsites;
    double d = 0.0;
    for (int i = 1; i < m.nsam; ++i)
        {
            d += 1. / static_cast<double>(i);
        }
    auto manual = S / d;
    BOOST_CHECK_CLOSE(w, manual, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_thetah)
// Simplest case
{
    auto h0 = Sequence::thetah(c, 0);
    auto h1 = Sequence::thetah(c, 1);
    auto m0 = manual_thetah(m, 0);
    auto m1 = manual_thetah(m, 1);
    BOOST_CHECK_CLOSE(h0, m0, 1e-6);
    BOOST_CHECK_CLOSE(h1, m1, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_faywuh)
// Fay and Wu's H is calculated using
// an aggregation distinct from
// thetah and thetapi.  That means
// we can compare results from the various functions
// in order to test.
{
    auto h = Sequence::thetah(c, 0);
    auto pi = Sequence::thetapi(c);
    auto fwh = Sequence::faywuh(c, 0);
    BOOST_CHECK_EQUAL(fwh + h, pi);
    BOOST_REQUIRE_EQUAL(pi - h, fwh);
}

BOOST_AUTO_TEST_CASE(test_thetah_multiple_derived_states)
{
    // Create a site with > 2 derived states
    auto f = Sequence::get_RowView(m, 0);
    for (std::size_t i = 0; i < f.size(); ++i)
        {
            f[i] = (i % 2 == 0.0) ? 2 : 3;
            if (i % 3 == 0.0)
                {
                    f[i] = 4;
                }
        }
    //We have to make data copies here so that 
    //max_allele is reset.
    Sequence::VariantMatrix m2(m.data, m.positions);
    BOOST_REQUIRE_THROW(auto h = Sequence::thetah(Sequence::AlleleCountMatrix(m2), 0), std::runtime_error);
    for (std::size_t i = 0; i < f.size(); ++i)
        {
            f[i] = 3;
        }
    Sequence::VariantMatrix m3(m.data, m.positions);
    BOOST_REQUIRE_NO_THROW(auto h = Sequence::thetah(Sequence::AlleleCountMatrix(m3), 0));
}

BOOST_AUTO_TEST_CASE(test_num_haplotypes)
{
    auto nh = Sequence::number_of_haplotypes(m);
    BOOST_REQUIRE_EQUAL(nh, 5);
}

BOOST_AUTO_TEST_CASE(test_unique_hap_at_any_index)
{
    for (std::size_t i = 0; i < m.nsam; ++i)
        {
            Sequence::VariantMatrix m2(m.data, m.positions);
            // We make a unique haplotype at this index in our copy of
            // the fixture
            auto cv = Sequence::get_ColView(m2, i);
            for (auto& j : cv)
                {
                    j = 2;
                }
            auto nh = Sequence::number_of_haplotypes(m2);
            std::vector<std::string> haps;
            for (std::size_t j = 0; j < m2.nsam; ++j)
                {
                    auto cvj = Sequence::get_ConstColView(m2, j);
                    std::string h;
                    for (auto state : cvj)
                        h += state;
                    haps.push_back(std::move(h));
                }
            std::set<std::string> uhaps(haps.begin(), haps.end());
            BOOST_REQUIRE_EQUAL(uhaps.size(), nh);
        }
}

BOOST_AUTO_TEST_CASE(test_haplotype_diversity)
{
    auto hd = Sequence::haplotype_diversity(m);

    // We calculate haplotype diversity here
    // in a brute-force manner:
    // 1. Build a vector of all haplotypes.
    // 2. Construct a vector of all unique haplotypes
    // 3. Explicitly cound the number of occurrences
    //    of each unique haplotype.
    std::vector<std::vector<std::int8_t>> haps;
    for (std::size_t i = 0; i < m.nsam; ++i)
        {
            auto col = Sequence::get_ConstColView(m, i);
            haps.emplace_back(
                std::vector<std::int8_t>(col.begin(), col.end()));
        }
    decltype(haps) uhaps;
    for (auto& h : haps)
        {
            if (std::find(uhaps.begin(), uhaps.end(), h) == uhaps.end())
                {
                    uhaps.push_back(h);
                }
        }
    double mhd = 0.0;
    for (auto& uh : uhaps)
        {
            auto c = std::count(haps.begin(), haps.end(), uh);
            mhd += static_cast<double>(c * (m.nsam - c));
        }
    mhd /= static_cast<double>(m.nsam * (m.nsam - 1));

    BOOST_CHECK_CLOSE(hd, mhd, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_rmin)
{
    auto rm = Sequence::rmin(m);
    BOOST_REQUIRE_EQUAL(rm, 1);
}

BOOST_AUTO_TEST_SUITE_END()
