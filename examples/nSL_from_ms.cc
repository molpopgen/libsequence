/*! \include nsl_from_ms.cc */
#include <cmath>
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <Sequence/summstats/nsl.hpp>

int
main(int argc, char** argv)
{
    auto vm = Sequence::from_msformat(std::cin);
    for (std::size_t i = 0; i < vm.nsites; ++i)
        {
            auto n = Sequence::nsl(vm, i, 0);
            if (!std::isnan(n.nsl))
                {
                    std::cout << vm.positions[i] << ' ' << n.nsl << ' '
                              << n.ihs << ' ' << n.core_count << '\n';
                }
        }
}
