//! \file VariantMatrixTest.cc @brief Tests for Sequence/VariantMatrix.hpp

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/windows.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <numeric> //for std::iota
#include <iterator>
#include "VariantMatrixFixture.hpp"

struct vmatrix_fixture
{
    std::vector<std::int8_t> input_data;
    std::vector<double> input_pos;
    Sequence::VariantMatrix m, m2;
    Sequence::AlleleCountMatrix c, c2;
    vmatrix_fixture()
        : input_data(make_input_data()), input_pos(make_intput_pos()),
          m(input_data, input_pos), m2(input_data, input_pos), c{ m }, c2{ m2 }
    {
        // The two VariantMatrix objects
        // have same data, but different internal
        // dimensions
        std::swap(m2.nsites, m2.nsam);
        m2.positions.resize(m2.nsites);
        std::iota(std::begin(m2.positions), std::end(m2.positions), 0.);
    }

    std::vector<std::int8_t>
    make_input_data()
    {
        int nsam = 20;
        int nsites = 5;
        std::vector<std::int8_t> rv;
        for (int i = 0; i < nsites; ++i)
            {
                for (int j = 0; j < nsam; ++j)
                    {
                        std::int8_t state = (j % 2 == 0.) ? 1 : 0;
                        rv.push_back(state);
                    }
            }
        return rv;
    }

    std::vector<double>
    make_intput_pos()
    {
        std::vector<double> rv;
        rv.resize(5);
        std::iota(rv.begin(), rv.end(), 0.);
        return rv;
    }
};

BOOST_FIXTURE_TEST_SUITE(VariantMatrixTest, vmatrix_fixture)

BOOST_AUTO_TEST_CASE(test_construction)
{
    Sequence::VariantMatrix m(std::vector<std::int8_t>(100, 1),
                              std::vector<double>(5, 0.0));
    BOOST_REQUIRE_EQUAL(m.nsites, 5);
    BOOST_REQUIRE_EQUAL(m.nsam, 20);
    BOOST_REQUIRE_EQUAL(m.max_allele, 1);
}

BOOST_AUTO_TEST_CASE(test_max_allele)
{
    m.data[3] = 5;
    Sequence::VariantMatrix vm(m.data, m.positions);
    BOOST_CHECK_EQUAL(vm.max_allele, 5);
    Sequence::AlleleCountMatrix vmc(vm);
    BOOST_REQUIRE_EQUAL(vmc.ncol, 6);
}

BOOST_AUTO_TEST_CASE(test_range_exceptions)
{
    BOOST_REQUIRE_THROW(m.at(m.nsites + 1, 0), std::out_of_range);
    BOOST_REQUIRE_THROW(m.at(0, m.nsam + 1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(test_iteration)
{
    for (std::size_t i = 0; i < m.nsam; ++i)
        {
            for (std::size_t j = 0; j < m.nsites; ++j)
                {
                    auto x = m.get(j, i);
                    std::int8_t ex = (i % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(x),
                                        static_cast<int>(ex));
                }
        }
}

BOOST_AUTO_TEST_CASE(test_bad_row_swap)
{
    auto a = Sequence::get_RowView(m, 0);
    auto b = Sequence::get_RowView(m2, 0);
    BOOST_REQUIRE_THROW(swap(a, b), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_bad_column_swap)
{
    auto a = Sequence::get_ColView(m, 0);
    auto b = Sequence::get_ColView(m2, 0);
    BOOST_REQUIRE_THROW(swap(a, b), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_row_views)
{
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto x = Sequence::get_RowView(m, i);
            for (auto j = x.begin(); j != x.end(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = x.cbegin(); j != x.cend(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.cbegin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = std::begin(x); j != std::end(x); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (std::size_t j = 0; j < x.size(); ++j)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(x[j]),
                                        static_cast<int>(ex));
                }
            std::size_t j = 0;
            for (auto xj : x)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(xj),
                                        static_cast<int>(ex));
                    ++j;
                }
        }
}

BOOST_AUTO_TEST_CASE(test_const_row_views)
{
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto x = Sequence::get_ConstRowView(m, i);

            for (auto j = x.begin(); j != x.end(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = x.cbegin(); j != x.cend(); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.cbegin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (auto j = std::begin(x); j != std::end(x); ++j)
                {
                    std::int8_t ex
                        = (std::distance(x.begin(), j) % 2 == 0.0) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(*j),
                                        static_cast<int>(ex));
                }
            for (std::size_t j = 0; j < x.size(); ++j)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(x[j]),
                                        static_cast<int>(ex));
                }
            std::size_t j = 0;
            for (auto xj : x)
                {
                    std::int8_t ex = (j % 2 == 0.) ? 1 : 0;
                    BOOST_REQUIRE_EQUAL(static_cast<int>(xj),
                                        static_cast<int>(ex));
                    ++j;
                }
        }
}

BOOST_AUTO_TEST_CASE(test_row_view_exceptions)
{
    BOOST_REQUIRE_THROW(Sequence::get_RowView(m, m.nsites + 1),
                        std::exception);
    BOOST_REQUIRE_THROW(Sequence::get_RowView(m, m.nsites + 1),
                        std::out_of_range);

    auto r = Sequence::get_RowView(m, 0);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::exception);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(test_const_row_view_exceptions)
{
    BOOST_REQUIRE_THROW(Sequence::get_ConstRowView(m, m.nsites + 1),
                        std::exception);
    BOOST_REQUIRE_THROW(Sequence::get_ConstRowView(m, m.nsites + 1),
                        std::out_of_range);

    auto r = Sequence::get_ConstRowView(m, 0);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::exception);
    BOOST_REQUIRE_THROW(r.at(m.nsam + 1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(test_row_view_iterators)
{
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto row = Sequence::get_RowView(m, i);
            BOOST_REQUIRE_EQUAL(std::distance(row.begin(), row.end()), m.nsam);
            BOOST_REQUIRE_EQUAL(std::distance(row.cbegin(), row.cend()),
                                m.nsam);
        }
}

BOOST_AUTO_TEST_CASE(test_column_views)
{
    for (std::size_t i = 0; i < m.nsam; ++i)
        {
            auto col = Sequence::get_ColView(m, i);
            std::int8_t state = (i % 2 == 0) ? 1 : 0;
            BOOST_REQUIRE_EQUAL(
                std::count(std::begin(col), std::end(col), !state), 0);
            BOOST_REQUIRE_EQUAL(std::count(col.rbegin(), col.rend(), !state),
                                0);
            BOOST_REQUIRE_EQUAL(std::count(col.cbegin(), col.cend(), !state),
                                0);
            BOOST_REQUIRE_EQUAL(std::count(col.crbegin(), col.crend(), !state),
                                0);

            BOOST_REQUIRE_EQUAL(std::distance(std::begin(col), std::end(col)),
                                m.nsites);
            BOOST_REQUIRE_EQUAL(std::distance(col.rbegin(), col.rend()),
                                m.nsites);

            // Check that iterators and reverse iterators have the expected
            // relationships:
            auto fwd = col.begin();
            auto rev = col.rbegin();
            for (; rev < col.rend(); ++rev)
                {
                    auto rf = std::distance(fwd, rev.base());
                    auto rb = std::distance(rev, col.rend());
                    BOOST_REQUIRE_EQUAL(rf, rb);
                }

            auto cfwd = col.cbegin();
            auto crev = col.crbegin();
            for (; crev < col.crend(); ++crev)
                {
                    auto rf = std::distance(cfwd, crev.base());
                    auto rb = std::distance(crev, col.crend());
                    BOOST_REQUIRE_EQUAL(rf, rb);
                }
        }
}

BOOST_AUTO_TEST_CASE(tesl_col_view_iterator_increment)
{
    auto x = Sequence::get_ConstColView(m, 0);
    auto b = x.begin();
    unsigned num_increments = 0;
    while (b < x.end())
        {
            b = b + 2;
            ++num_increments;
        }
    BOOST_REQUIRE_EQUAL(num_increments, 3);
}

BOOST_AUTO_TEST_CASE(test_column_view_invalid_compare)
{
    auto c0 = Sequence::get_ConstColView(m, 0);
    auto c1 = Sequence::get_ConstColView(m, 1);
    BOOST_REQUIRE_NO_THROW(std::distance(c0.begin(), c0.end()));
    BOOST_REQUIRE_THROW(std::distance(c0.begin(), c1.begin()),
                        std::invalid_argument);
}

// The remaining tests apply STL algorithms to column iterators,
// which is a good stress test.  We've already done count above.

BOOST_AUTO_TEST_CASE(test_accumulate)
{
    auto c = Sequence::get_ConstColView(m, 0);
    int sum = static_cast<int>(std::accumulate(c.cbegin(), c.cend(), 0));
    BOOST_REQUIRE_EQUAL(sum, static_cast<int>(m.nsites));
    c = Sequence::get_ConstColView(m, 1);
    sum = static_cast<int>(std::accumulate(c.cbegin(), c.cend(), 0));
    BOOST_REQUIRE_EQUAL(sum, 0);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(VariantMatrixWindowTest, dataset)

BOOST_AUTO_TEST_CASE(tests_windows_size_0)
{
	auto w = Sequence::make_window(m, 10, 10);
}

BOOST_AUTO_TEST_CASE(tests_windows_size_1)
{
    for (std::size_t i = 0; i < m.positions.size(); ++i)
        {
            auto w = Sequence::make_window(m, m.positions[i], m.positions[i]);
            BOOST_REQUIRE_EQUAL(w.positions[0], m.positions[i]);
            auto r = Sequence::get_ConstRowView(m, i);
            BOOST_REQUIRE_EQUAL(
                std::mismatch(w.data.begin(), w.data.end(), r.begin()).first
                    == w.data.end(),
                true);
        }
}

BOOST_AUTO_TEST_CASE(test_windows_multi_site)
{
	auto w = Sequence::make_window(m, 0.11, 0.29);
	BOOST_REQUIRE_EQUAL(w.nsites, 1);
	BOOST_REQUIRE_EQUAL(w.positions[0], 0.2);
}

BOOST_AUTO_TEST_SUITE_END()
