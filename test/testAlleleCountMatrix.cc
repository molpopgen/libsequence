//! \file testAlleleCountMatrix.cc @brief Tests for Sequence/VariantMatrix.hpp
#include "VariantMatrixFixture.hpp"
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <numeric> //for std::iota
#include <iterator>

BOOST_FIXTURE_TEST_SUITE(test_allele_count_matrix, dataset)

BOOST_AUTO_TEST_CASE(test_max_allele_exception)
{
	//Change some data in m so that m[i] > m.max_allele
	m.data[0] = 5;

	BOOST_REQUIRE_THROW(Sequence::AlleleCountMatrix ac(m),std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

