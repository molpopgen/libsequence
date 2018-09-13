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

} // namespace Sequence

#endif
