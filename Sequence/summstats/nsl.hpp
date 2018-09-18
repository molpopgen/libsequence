/// \file Sequence/summstats/nsl.hpp
/// \brief nSL and iHS
#ifndef SEQUENCE_SUMMSTATS_NSL_HPP__
#define SEQUENCE_SUMMSTATS_NSL_HPP__

#include <vector>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct nSLiHS
    /// Stores the results of nSL and iHS calculations.
    /// See Sequence::nsl for details.
    /// \ingroup popgenanalysis
    {
        /// The nSL statistic \cite Ferrer-Admetlla2014-wa
        double nsl;
        /// The iHS statistic, calculated according to \cite Ferrer-Admetlla2014-wa
        double ihs;
        /// Count of non-reference,
        /// non-missing allele.
        std::int32_t core_count;
    };

    /*! \brief nSL and iHS statistics
     * \param m A VariantMatrix
     * \param core The index of the core site
     * \param refstate The value of the reference/ancestral allelic state
     *
     * \return an nSLiHS object
     * \ingroup popgenanalysis
     *
     * See nSL_from_ms.cc for example
     *
     * See \cite Ferrer-Admetlla2014-wa for details.
     */
    nSLiHS nsl(const VariantMatrix& m, const std::size_t core,
               const std::int8_t refstate);

    /*! \brief nSL and iHS statistics
     * \param m A VariantMatrix
     * \param refstate The value of the reference/ancestral allelic state
     *
     * \return vector of nSLiHS objects (one for each site)
     * \ingroup popgenanalysis
     *
     * This function differs from the version working 
     * on a core site in that it uses an efficient matrix-based
     * method to dynamically update suffix lengths as each 
     * core site is processed.  The result is a huge runtime
     * reduction compared to calculating the statistic
     * for each core site on its own.
     *
     * See \cite Ferrer-Admetlla2014-wa for details.
     */
        std::vector<nSLiHS> nsl(const VariantMatrix& m,
                            const std::int8_t refstate);
} // namespace Sequence

#endif
