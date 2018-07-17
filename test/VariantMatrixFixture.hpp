#ifndef LIBSEQUENCE_TESTS_VARIANTMATRIXFIXTURE_HPP__
#define LIBSEQUENCE_TESTS_VARIANTMATRIXFIXTURE_HPP__

#include <Sequence/VariantMatrix.hpp>
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

#endif
