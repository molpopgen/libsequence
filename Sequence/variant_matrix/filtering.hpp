#ifndef SEQUENCE_VARIANT_MATRIX_FILTERING_HPP_
#define SEQUENCE_VARIANT_MATRIX_FILTERING_HPP_

#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <functional>
#include <cstdint>

namespace Sequence
{
    std::int32_t filter_sites(VariantMatrix &m,
                              const std::function<bool(const RowView &)> &f);

    std::int32_t
    filter_haplotypes(VariantMatrix &m,
                      const std::function<bool(const ColView &)> &f);

    std::int32_t filter_sites(VariantMatrix &m,
                              const std::function<bool(const ConstRowView &)> &f);

    std::int32_t
    filter_haplotypes(VariantMatrix &m,
                      const std::function<bool(const ConstColView &)> &f);
}

#endif
