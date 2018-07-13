//! \file testClassicSummstats.cc @brief unit tests for Sequence::ClustalW and Sequence::phylipData

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
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
                        1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
             positions_type{ 0.1, 0.2, 0.3 } }
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(test_classic_stats, dataset)

BOOST_AUTO_TEST_CASE(test_thetapi)
{
    auto pi = Sequence::thetapi(m);

    double manual = 0.0;
    for (std::size_t i = 0; i < m.nsites; ++i)
        {
            auto r = Sequence::get_ConstRowView(m, i);
            int ndiffs = 0, ncomps = 0;
            for (std::size_t j = 0; j < m.nsam - 1; ++j)
                {
                    for (std::size_t k = j + 1; k < m.nsam ; ++k)
                        {
                            if (r[j] != r[k])
                                ++ndiffs;
                            ++ncomps;
                        }
                }
            manual
                += static_cast<double>(ndiffs) / static_cast<double>(ncomps);
        }
    // Cannot require equal b/c we aren't doing ops
    // in same order.
    BOOST_CHECK_CLOSE(pi, manual, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
