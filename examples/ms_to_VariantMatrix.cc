/*! \include ms_to_VariantMatrix.cc */
#include <iostream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>

int
main(int argc, char** argv)
{
    do
        {
            auto vm = Sequence::from_msformat(std::cin);
            Sequence::to_msformat(vm, std::cout);
            std::cout << '\n';
        }
    while (!std::cin.eof());
}
