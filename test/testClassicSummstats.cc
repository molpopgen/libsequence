//! \file testClassicSummstats.cc @brief unit tests for Sequence::ClustalW and Sequence::phylipData

#include <cmath>
#include <algorithm>
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/classics.hpp>
#include <boost/test/unit_test.hpp>

struct dataset
{
    using data_type = decltype(Sequence::VariantMatrix::data);
    using positions_type = decltype(Sequence::VariantMatrix::positions);
    Sequence::VariantMatrix m;
    dataset()
        : m{ data_type{ 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
                        1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
             positions_type{ 0.1, 0.2, 0.3 } }
    {
    }
};

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
    auto pi = Sequence::thetapi(m);

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
    auto pi = Sequence::thetapi(m);

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

BOOST_AUTO_TEST_CASE(test_thetaw)
// The above two tests imply that
// thetaw is correct.  Thus,
// all we're doing below is
// checking the denominator.
{
    auto w = Sequence::thetaw(m);
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
    auto h0 = Sequence::thetah(m, 0);
    auto h1 = Sequence::thetah(m, 1);
    auto m0 = manual_thetah(m, 0);
    auto m1 = manual_thetah(m, 1);
    BOOST_CHECK_CLOSE(h0, m0, 1e-6);
    BOOST_CHECK_CLOSE(h1, m1, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
