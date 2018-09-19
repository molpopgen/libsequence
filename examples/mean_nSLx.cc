/*! \include mean_nSLx.cc */
#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <Sequence/summstats/nslx.hpp>

int
main(int argc, char** argv)
{
    int x = std::atoi(argv[1]);
    while (!std::cin.eof())
        {
            auto vm = Sequence::from_msformat(std::cin);
            auto nsl_stats = Sequence::nslx(vm, 0, x);
            double sum = 0.0;
            unsigned n = 0;
            for (auto& i : nsl_stats)
                {
                    if (std::isfinite(i.nsl))
                        {
                            sum += i.nsl;
                            ++n;
                        }
                }
            std::cout << sum / static_cast<double>(n) << '\n';
        }
}
