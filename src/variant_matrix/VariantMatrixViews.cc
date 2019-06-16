#include <Sequence/VariantMatrixViews.hpp>
#include <stdexcept>

namespace
{
    template <typename T, typename VM>
    T
    row_view_wrapper(VM& m, const std::size_t row)
    {
        if (row >= m.nsites())
            {
                throw std::out_of_range("row index out of range");
            }
        auto x = m.data() + m.genotype_row_offset() * m.nsam()
                 + m.genotype_col_offset() + row * m.nsam();
        return T(x, m.nsam());
    }

    template <typename T, typename VM>
    T
    col_view_wrapper(VM& m, const std::size_t col)
    {
        if (col >= m.nsam())
            {
                throw std::out_of_range("column index out of range");
            }
        auto x = m.data() + m.genotype_row_offset() * m.genotype_stride()
                 + m.genotype_col_offset() + col;
        return T(x, m.genotype_stride() * m.nsites(), m.genotype_stride());
    }

    template <typename T, typename VM>
    T
    const_row_view_wrapper(VM& m, const std::size_t row)
    {
        if (row >= m.nsites())
            {
                throw std::out_of_range("row index out of range");
            }
        auto x = m.cdata() + m.genotype_row_offset() * m.nsam()
                 + m.genotype_col_offset() + row * m.nsam();
        return T(x, m.nsam());
    }

    template <typename T, typename VM>
    T
    const_col_view_wrapper(VM& m, const std::size_t col)
    {
        if (col >= m.nsam())
            {
                throw std::out_of_range("column index out of range");
            }
        auto x = m.cdata() + m.genotype_row_offset() * m.genotype_stride()
                 + m.genotype_col_offset() + col;
        return T(x, m.genotype_stride() * m.nsites(), m.genotype_stride());
    }
} // namespace

namespace Sequence
{

    ConstRowView
    get_RowView(const VariantMatrix& m, const std::size_t row)
    {
        return row_view_wrapper<ConstRowView>(m, row);
    }

    RowView
    get_RowView(VariantMatrix& m, const std::size_t row)
    {
        return row_view_wrapper<RowView>(m, row);
    }

    ConstRowView
    get_ConstRowView(const VariantMatrix& m, const std::size_t row)
    {
        return const_row_view_wrapper<ConstRowView>(m, row);
    }

    ConstRowView
    get_ConstRowView(VariantMatrix& m, const std::size_t row)
    {
        return const_row_view_wrapper<ConstRowView>(m, row);
    }

    ColView
    get_ColView(VariantMatrix& m, const std::size_t col)
    {
        return col_view_wrapper<ColView>(m, col);
    }
    ConstColView

    get_ColView(const VariantMatrix& m, const std::size_t col)
    {
        return col_view_wrapper<ConstColView>(m, col);
    }

    ConstColView
    get_ConstColView(VariantMatrix& m, const std::size_t col)
    {
        return const_col_view_wrapper<ConstColView>(m, col);
    }

    ConstColView
    get_ConstColView(const VariantMatrix& m, const std::size_t col)
    {
        return const_col_view_wrapper<ConstColView>(m, col);
    }
} // namespace Sequence
