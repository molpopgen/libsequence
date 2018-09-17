/*! \include nSL_from_ms.cc */
#include <cmath>
#include <iostream>
#include <cassert>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <Sequence/summstats/nsl.hpp>
#include <Sequence/summstats/nslx.hpp>

int
main(int argc, char** argv)
{
    int x = std::atoi(argv[1]);
    auto vm = Sequence::from_msformat(std::cin);
    auto nsl_stats = Sequence::nslx(vm, 0, x);
    for(auto & s : nsl_stats){std::cout << s.nsl << ' ' << s.ihs << ' ' << s.core_count << '\n'; }
}

