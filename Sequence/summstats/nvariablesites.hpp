/// \file Sequence/summstats/nvariablesites.hpp
/// \brief Calculate total numbers of polymorphisms
#ifndef SEQUENCE_SUMMSTATS_NVARIABLESITES_HPP__
#define SEQUENCE_SUMMSTATS_NVARIABLESITES_HPP__

#include <cstdint>
#include <Sequence/AlleleCountMatrix.hpp>

namespace Sequence
{
    /*! \brief Number of polymorphic sites
     *
     * Returns the number of sites with more than one non-missing state
     * \param m An AlleleCountMatrix
     * \return std::uint32_t
     * \ingroup popgenanalysis
     */
    std::uint32_t nvariable_sites(const AlleleCountMatrix& m);

    /*! \brief Number of bi-allelic sites
     *
     * Return the number of sites with exactly two non-missing states.
     * \param m An AlleleCountMatrix
     * \return std::uint32_t
     * \ingroup popgenanalysis
     */
    std::uint32_t nbiallelic_sites(const AlleleCountMatrix& m);

    /*! \brief Total number of mutations in the sample
     *
     * Return \f$\sum_{i=0}^{i=m.nsites-1}I(i)\f$ where \f$I(i)\f$
     * is \f$k_i - 1\f$ if \f$k_i\f$, the number of states at the \f$i^{th}\f$ site,
     * is greater than one, and zero otherwise.
     *
     * \param m An AlleleCountMatrix
     * \return std::uint32_t
     * \ingroup popgenanalysis
     */
    std::uint32_t total_number_of_mutations(const AlleleCountMatrix& m);
} // namespace Sequence

#endif
