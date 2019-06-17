//! \file VariantMatrixTest.cc @brief Tests for Sequence/VariantMatrix.hpp

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/windows.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <numeric> //for std::iota
#include <iterator>
#include <sstream>
#include <fstream>
#include <iostream>
#include "msprime_data_fixture.hpp"

bool
is_const_row_view(Sequence::ConstRowView& v)
{
    return true;
}

bool
is_const_row_view(Sequence::RowView& v)
{
    return false;
}

bool
is_const_col_view(Sequence::ConstColView& v)
{
    return true;
}

bool
is_const_col_view(Sequence::ColView& v)
{
    return false;
}

BOOST_AUTO_TEST_SUITE(BasicVariantMatrixTests)

BOOST_AUTO_TEST_CASE(test_construction)
{
    Sequence::VariantMatrix m(std::vector<std::int8_t>(100, 1),
                              std::vector<double>(5, 0.0));
    BOOST_REQUIRE_EQUAL(m.nsites(), 5);
    BOOST_REQUIRE_EQUAL(m.nsam(), 20);
    BOOST_REQUIRE_EQUAL(m.max_allele(), 1);
}

BOOST_AUTO_TEST_CASE(construct_empty_VariantMatrix_from_move)
{
    std::vector<std::int8_t> d{};
    std::vector<double> p{};
    Sequence::VariantMatrix m(std::move(d), std::move(p));
    BOOST_REQUIRE_EQUAL(m.nsam(), 0);
    BOOST_REQUIRE_EQUAL(m.nsites(), 0);
}

BOOST_AUTO_TEST_CASE(construct_empty_VariantMatrix_from_init_lists)
{
    Sequence::VariantMatrix m(std::vector<std::int8_t>{},
                              std::vector<double>{});
    BOOST_REQUIRE_EQUAL(m.nsam(), 0);
    BOOST_REQUIRE_EQUAL(m.nsites(), 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(VariantMatrixTest, vmatrix_from_msprime)


BOOST_AUTO_TEST_CASE(test_max_allele)
{
    m.data()[3] = 5;
    Sequence::VariantMatrix vm(
        std::vector<std::int8_t>(m.data(), m.data() + m.nsites() * m.nsam()),
        std::vector<double>(m.pbegin(), m.pend()));
    BOOST_CHECK_EQUAL(vm.max_allele(), 5);
    Sequence::AlleleCountMatrix vmc(vm);
    BOOST_REQUIRE_EQUAL(vmc.ncol, 6);
}

BOOST_AUTO_TEST_CASE(test_range_exceptions)
{
    BOOST_REQUIRE_THROW(m.at(m.nsites() + 1, 0), std::out_of_range);
    BOOST_REQUIRE_THROW(m.at(0, m.nsam() + 1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(test_row_views)
{
    for (std::size_t i = 0; i < m.nsites(); ++i)
        {
            auto x = Sequence::get_RowView(m, i);
            BOOST_REQUIRE_EQUAL(is_const_row_view(x), false);
        }
}

BOOST_AUTO_TEST_CASE(test_const_row_views)
{
    for (std::size_t i = 0; i < m.nsites(); ++i)
        {
            auto x = Sequence::get_ConstRowView(m, i);
            BOOST_REQUIRE_EQUAL(is_const_row_view(x), true);
        }
}

BOOST_AUTO_TEST_CASE(test_row_view_exceptions)
{
    BOOST_REQUIRE_THROW(Sequence::get_RowView(m, m.nsites() + 1),
                        std::exception);
    BOOST_REQUIRE_THROW(Sequence::get_RowView(m, m.nsites() + 1),
                        std::out_of_range);

    auto r = Sequence::get_RowView(m, 0);
    BOOST_REQUIRE_THROW(r.at(m.nsam() + 1), std::exception);
    BOOST_REQUIRE_THROW(r.at(m.nsam() + 1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(test_const_row_view_exceptions)
{
    BOOST_REQUIRE_THROW(Sequence::get_ConstRowView(m, m.nsites() + 1),
                        std::exception);
    BOOST_REQUIRE_THROW(Sequence::get_ConstRowView(m, m.nsites() + 1),
                        std::out_of_range);

    auto r = Sequence::get_ConstRowView(m, 0);
    BOOST_REQUIRE_THROW(r.at(m.nsam() + 1), std::exception);
    BOOST_REQUIRE_THROW(r.at(m.nsam() + 1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(test_row_view_iterators)
{
    for (std::size_t i = 0; i < m.nsites(); ++i)
        {
            auto row = Sequence::get_RowView(m, i);
            BOOST_REQUIRE_EQUAL(std::distance(row.begin(), row.end()),
                                m.nsam());
            BOOST_REQUIRE_EQUAL(std::distance(row.cbegin(), row.cend()),
                                m.nsam());
        }
}

BOOST_AUTO_TEST_CASE(test_column_views)
{
    for (std::size_t i = 0; i < m.nsam(); ++i)
        {
            auto col = Sequence::get_ColView(m, i);
            BOOST_REQUIRE_EQUAL(is_const_col_view(col), false);
        }
}

BOOST_AUTO_TEST_CASE(test_column_view_invalid_compare)
{
    auto c0 = Sequence::get_ConstColView(m, 0);
    auto c1 = Sequence::get_ConstColView(m, 1);
    BOOST_REQUIRE_NO_THROW(std::distance(c0.begin(), c0.end()));
    BOOST_REQUIRE_THROW(std::distance(c0.begin(), c1.begin()),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_deep_copy, vmatrix_from_msprime)

BOOST_AUTO_TEST_CASE(test_variant_matrix_deepcopy)
{
    auto c = m.deepcopy();
    BOOST_REQUIRE_EQUAL(m == c, true);
    BOOST_REQUIRE_EQUAL(m != c, false);
}

BOOST_AUTO_TEST_CASE(test_row_swap)
{
    auto m2 = m.deepcopy();
    auto a = Sequence::get_RowView(m, 0);
    auto b = Sequence::get_RowView(m2, 0);
    BOOST_REQUIRE_NO_THROW(swap(a, b));
}

BOOST_AUTO_TEST_CASE(test_column_swap)
{
    auto m2 = m.deepcopy();
    auto a = Sequence::get_ColView(m, 0);
    auto b = Sequence::get_ColView(m2, 0);
    BOOST_REQUIRE_NO_THROW(swap(a, b));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_msformat_round_trips, vmatrix_from_msprime)

BOOST_AUTO_TEST_CASE(test_io_roundtrip)
{
    std::ostringstream o;
    Sequence::to_msformat(m, o);
    std::istringstream in(o.str());
    auto vm = Sequence::from_msformat(in);
    BOOST_REQUIRE_EQUAL(m == vm, true);
    BOOST_REQUIRE_EQUAL(m != vm, false);
}

BOOST_AUTO_TEST_SUITE_END()
