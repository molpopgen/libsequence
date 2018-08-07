/// \file Sequence/summstats/thetal.hpp
/// \brief Zeng et al. \f$\hat\theta_L\f$
#ifndef SEQUENCE_SUMMSTATS_THETAL_HPP__
#define SEQUENCE_SUMMSTATS_THETAL_HPP__

#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    /*! \brief Zeng et al. \f$\hat\theta_L\f$
     * \param m a VariantMatrix
     * \param refstate The ancestral state
     * \return double
     *
     * See \cite Zeng2006-is for details.
     * \ingroup popgenanalysis
     */
    double thetal(const VariantMatrix& m, const std::int8_t refstate);

    /*! \brief Zeng et al. \f$\hat\theta_L\f$
     * \param m a VariantMatrix
     * \param refstate Vector of ancestral states.
     * \return double
     *
     * See \cite Zeng2006-is for details.
     * \ingroup popgenanalysis
     */
    double thetal(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);
    double thetal(const AlleleCountMatrix& ac, const std::int8_t refstate);
    double thetal(const AlleleCountMatrix& m,
                  const std::vector<std::int8_t>& refstates);
} // namespace Sequence

#endif
