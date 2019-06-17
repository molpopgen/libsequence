#ifndef LIBSEQUENCE_TESTS_VARIANTMATRIXFIXTURE_HPP
#define LIBSEQUENCE_TESTS_VARIANTMATRIXFIXTURE_HPP

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/AlleleCountMatrix.hpp>

struct invariantdataset
{
    using data_type = std::vector<std::int8_t>;
    using positions_type = std::vector<double>;
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
