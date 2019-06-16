#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/windows.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <sstream>
#include <iostream>
#include "msformatdata.hpp"

BOOST_AUTO_TEST_SUITE(testVariantMatrixWindows)

BOOST_AUTO_TEST_CASE(test_windows)
{
    std::istringstream i(get_msformat_data());
    auto vm = Sequence::from_msformat(i);
    Sequence::to_msformat(vm, std::cout);
}

BOOST_AUTO_TEST_SUITE_END()
