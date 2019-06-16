/*! \include ms_to_VariantMatrix.cc */
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <Sequence/variant_matrix/windows.hpp>

int
main(int argc, char** argv)
{
    do
        {
            auto vm = Sequence::from_msformat(std::cin);
            Sequence::to_msformat(vm, std::cout);
            std::cout << '\n';
            for (double i = 0; i < 1. - 0.01; i += 0.1)
                {
                    auto w = Sequence::make_window(vm, i, i + 0.1);
                    Sequence::to_msformat(w, std::cout);
                    std::cout << '\n';
                    w = Sequence::make_slice(vm, i, i + 0.1, 3, 8);
                    std::cout << "-------\n";
                    for (std::size_t k = 0; k < w.nsites(); ++k)
                        {
                            auto c = Sequence::get_ConstRowView(w, k);
                            for (auto l = 0; l < c.size(); ++l)
                                {
                                    std::cout << int(c[l]);
                                }
                            std::cout << '\n';
                        }
                    std::cout << "-------\n";
                    std::cerr << "gst = " << w.genotype_stride() << ' '
                              << w.nsites() << ' ' << w.nsam() << '\n';
                    Sequence::to_msformat(w, std::cout);
                    std::cout << '\n';
                }
        }
    while (!std::cin.eof());
}
