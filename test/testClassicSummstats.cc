//! \file testClassicSummstats.cc @brief unit tests for Sequence::ClustalW and Sequence::phylipData

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/summstats/classics.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <unistd.h>

struct dataset
{
    using data_type = decltype(Sequence::VariantMatrix::data);
    using positions_type = decltype(Sequence::VariantMatrix::positions);
    Sequence::VariantMatrix m;
    dataset()
        : m{ data_type{ 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
               1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1 },
             positions_type{ 0.1, 0.2, 0.3 } }
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(test_classic_stats, dataset)

BOOST_AUTO_TEST_CASE(test_thetapi)
{
    auto pi = Sequence::thetapi(m);
}

BOOST_AUTO_TEST_SUITE_END()
