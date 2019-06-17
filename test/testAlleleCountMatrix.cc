//! \file testAlleleCountMatrix.cc @brief Tests for Sequence/VariantMatrix.hpp
#include "msprime_data_fixture.hpp"
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/variant_matrix/windows.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <numeric> //for std::iota
#include <iterator>

BOOST_FIXTURE_TEST_SUITE(test_allele_count_matrix, vmatrix_from_msprime)

BOOST_AUTO_TEST_CASE(test_max_allele_exception)
{
    //Change some data in m so that m[i] > m.max_allele
    m.data()[0] = 5;

    BOOST_REQUIRE_THROW(Sequence::AlleleCountMatrix ac(m), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(counts_from_windows)
{
    for (std::size_t i = 0; i < m.nsites(); ++i)
        {
            auto w = Sequence::make_window(m, m.position(i), m.position(i));
            BOOST_REQUIRE_NO_THROW(Sequence::AlleleCountMatrix ac(w));
        }
}

BOOST_AUTO_TEST_SUITE_END()

