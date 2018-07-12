#ifndef SEQUENCE_VARIANT_MATRIX_VIEWS_HPP__
#define SEQUENCE_VARIANT_MATRIX_VIEWS_HPP__

#include "VariantMatrix.hpp"
#include "bits/variant_matrix_views_internal.hpp"

namespace Sequence
{
    /// \brief View of a VariantMatrix row (a variable site)
    /// \ingroup variantmatrix
    using RowView = internal::row_view_<std::int8_t*>;
    /// \brief Const view of a VariantMatrix row (a variable site)
    /// \ingroup variantmatrix
    using ConstRowView = internal::row_view_<const std::int8_t*>;
    /// \brief View of a VariantMatrix column ("haplotype")
    /// \ingroup variantmatrix
    using ColView = internal::col_view_<std::int8_t*>;
    /// \brief Const view of a VariantMatrix column ("haplotype")
    /// \ingroup variantmatrix
    using ConstColView = internal::col_view_<const std::int8_t*>;

    // Rather than have member functions, we will have standalone functions:

    // The following could (should?) be declared noexcept(false):

    /// \brief Return a ConstRowView from VariantMatrix `m` at row `row`.
    /// std::out_of_range is thrown if `row` is out of range.
    /// \ingroup variantmatrix
    ConstRowView get_RowView(const VariantMatrix& m, const std::size_t row);

    /// Return a RowView from VariantMatrix `m` at row `row`.
    /// std::out_of_range is thrown if `row` is out of range.
    /// \ingroup variantmatrix
    RowView get_RowView(VariantMatrix& m, const std::size_t row);

    /// \brief Return a ConstRowView from VariantMatrix `m` at row `row`.
    /// std::out_of_range is thrown if `row` is out of range.
    /// \ingroup variantmatrix
    ConstRowView get_ConstRowView(const VariantMatrix& m,
                                  const std::size_t row);

    /// \brief Return a ConstRowView from VariantMatrix `m` at row `row`.
    /// std::out_of_range is thrown if `row` is out of range.
    /// \ingroup variantmatrix
    ConstRowView get_ConstRowView(VariantMatrix& m, const std::size_t row);

    /// \brief Return a ColView from VariantMatrix `m` at col `col`.
    /// std::out_of_range is thcoln if `col` is out of range.
    /// \ingroup variantmatrix
    ColView get_ColView(VariantMatrix& m, const std::size_t col);

    /// \brief Return a ConstColView from VariantMatrix `m` at col `col`.
    /// std::out_of_range is thcoln if `col` is out of range.
    /// \ingroup variantmatrix
    ConstColView get_ColView(const VariantMatrix& m, const std::size_t col);

    /// \brief Return a ConstColView from VariantMatrix `m` at col `col`.
    /// std::out_of_range is thcoln if `col` is out of range.
    /// \ingroup variantmatrix
    ConstColView get_ConstColView(VariantMatrix& m, const std::size_t col);
    
    /// \brief Return a ConstColView from VariantMatrix `m` at col `col`.
    /// std::out_of_range is thcoln if `col` is out of range.
    /// \ingroup variantmatrix
    ConstColView get_ConstColView(const VariantMatrix& m,
                                  const std::size_t col);
}

#endif
