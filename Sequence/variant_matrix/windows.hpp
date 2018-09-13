#ifndef SEQUENCE_VARIANT_MATRIX_WINDOWS_HPP
#define SEQUENCE_VARIANT_MATRIX_WINDOWS_HPP

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    /*! \brief Return a window from a VariantMatrix
     * \param m A VariantMatrix
     * \param beg Beginning of window
     * \param end End of window
     *
     * \note The window intervals are open, [beg,end]
     */
    VariantMatrix make_window(const VariantMatrix& m, const double beg,
                              const double end);
    /*! \brief Return a slice from a VariantMatrix
     * \param m A VariantMatrix
     * \param beg Beginning of window
     * \param end End of window
     * \param i index of first haplotype to include
     * \param j one past last haplotype to include
     *
     * The result is a variant matrix including positions [beg,end]
     * and samples [i,j) from \a m.  Note that the sample interval is 
     * half-open!
     */
    
    VariantMatrix make_slice(const VariantMatrix& m, const double beg,
                             const double end,
                             const std::size_t i,
                             const std::size_t j);
} // namespace Sequence

#endif
