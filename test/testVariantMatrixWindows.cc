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
    for (double i = 0.0; i < 1.0 - 1e-4; i += 0.1)
        {
            auto w = Sequence::make_window(vm, i, i + 0.1);
            auto pb = std::lower_bound(vm.pbegin(), vm.pend(), i);
            std::size_t offset = pb - vm.pbegin();
            for (std::size_t site = 0; site < w.nsites(); ++site)
                {
                    auto window_site = Sequence::get_ConstRowView(w, site);
                    auto matrix_site
                        = Sequence::get_ConstRowView(vm, offset + site);
                    auto m
                        = std::mismatch(window_site.begin(), window_site.end(),
                                        matrix_site.begin());
                    BOOST_REQUIRE_EQUAL(m.first == window_site.end(), true);
                }
        }
}

BOOST_AUTO_TEST_SUITE_END()
