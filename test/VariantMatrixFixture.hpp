#ifndef LIBSEQUENCE_TESTS_VARIANTMATRIXFIXTURE_HPP
#define LIBSEQUENCE_TESTS_VARIANTMATRIXFIXTURE_HPP

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/AlleleCountMatrix.hpp>

struct dataset
{
    using data_type = decltype(Sequence::VariantMatrix::data);
    using positions_type = decltype(Sequence::VariantMatrix::positions);
    Sequence::VariantMatrix m;
    Sequence::AlleleCountMatrix c;
    dataset()
        : m{ data_type{ 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
                        1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
             positions_type{ 0.1, 0.2, 0.3 } },
          c{ m }
    {
    }
};

struct invariantdataset
{
    using data_type = decltype(Sequence::VariantMatrix::data);
    using positions_type = decltype(Sequence::VariantMatrix::positions);
    Sequence::VariantMatrix empty, invariant;
    Sequence::AlleleCountMatrix empty_counts, invariant_counts;
    invariantdataset()
        : empty{ data_type{}, positions_type{} },
          invariant{ data_type{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                     positions_type{ 0.1, 0.2, 0.3 } },
          empty_counts(empty), invariant_counts(invariant)
    {
    }
};

#endif
