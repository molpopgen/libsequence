/// \file Sequence/summstats/thetah.hpp
/// \brief Fay and Wu's \f$\hat\theta_H\f$.
#ifndef SEQUENCE_SUMMSTATS_THETAH_HPP__
#define SEQUENCE_SUMMSTATS_THETAH_HPP__

#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    /*! \brief Fay and Wu's \f$\hat\theta_H\f$.
     * \param m a VariantMatrix
     * \param refstate The ancestral state
     * \return double
     *
     * See \cite Fay2000-ef for details.
     * \ingroup popgenanalysis
     */
    double thetah(const VariantMatrix& m, const std::int8_t refstate);

    /*! \brief Fay and Wu's \f$\hat\theta_H\f$.
     * \param m a VariantMatrix
     * \param refstate Vector of ancestral states.
     * \return double
     *
     * See \cite Fay2000-ef for details.
     * \ingroup popgenanalysis
     */
    double thetah(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);
} // namespace Sequence

#endif
