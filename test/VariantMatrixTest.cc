//! \file ComparisonsTests.cc @brief Tests for Sequence/Comparisons.hpp
#define BOOST_TEST_MODULE VariantMatrixTest

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <iterator>
#include <iostream>
struct vmatrix_fixture
{
    std::vector<std::int8_t> input_data;
    std::vector<double> input_pos;
    Sequence::VariantMatrix m;

    vmatrix_fixture()
        : input_data(make_input_data()), input_pos(make_intput_pos()),
          m(input_data, input_pos)
    {
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

BOOST_FIXTURE_TEST_CASE(test_construction, vmatrix_fixture)
{
    Sequence::VariantMatrix m(std::vector<std::int8_t>(100, 1),
                              std::vector<double>(5, 0.0));
    BOOST_REQUIRE_EQUAL(m.nsites, 5);
    BOOST_REQUIRE_EQUAL(m.nsam, 20);
}

BOOST_FIXTURE_TEST_CASE(test_iteration, vmatrix_fixture)
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
