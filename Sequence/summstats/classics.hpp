/// \file Sequence/summstats/classics.hpp
/// \brief "Classic" summaries of variation data.
#ifndef SEQUENCE_SUMMSTATS_CLASSICS_HPP__
#define SEQUENCE_SUMMSTATS_CLASSICS_HPP__

#include <Sequence/VariantMatrix.hpp>

#include "thetapi.hpp"
#include "thetaw.hpp"
#include "thetah.hpp"
#include "thetal.hpp"
#include "nvariablesites.hpp"
#include "allele_counts.hpp"

namespace Sequence
{
    /*! \brief Tajima's D
     * \param m A VariantMatrix
     * \return Tajima's D
     * \note nan is returned if m is empty/invariant.
     *
     * See \cite Tajima1989-de for details.
     *
     * Included via Sequence/summstats.hpp or
     * Sequence/summstats/classics.hpp
     *
     * \ingroup popgenanalysis
     */
    double tajd(const VariantMatrix& m);

    double
    tajd(const AlleleCountMatrix& ac);

    /*! The H' statistic
     * \param m A VariantMatrix
     * \param refstate How the ancestral state is encoded.
     * \return H'
     * \note nan is returned if m is empty/invariant.  It is
     * an error to consider \a refstate as anything other than
     * the ancestral state.
     *
     * See \cite Zeng2006-is for details.
     *
     * Included via Sequence/summstats.hpp or
     * Sequence/summstats/classics.hpp
     *
     * \ingroup popgenanalysis
     */
    double hprime(const VariantMatrix& m, const std::int8_t refstate);

    /*! The H' statistic
     * \param m A VariantMatrix
     * \param refstates A vector of ancestral states, equal in length to m.sites
     * \return H'
     * \note nan is returned if m is empty/invariant.  It is
     * an error to consider \a refstates as anything other than
     * the ancestral state at each site.
     *
     * See \cite Zeng2006-is for details.
     *
     * Included via Sequence/summstats.hpp or
     * Sequence/summstats/classics.hpp
     *
     * \ingroup popgenanalysis
     */
    double hprime(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);

    /*! \brief Fay and Wu's H.
     * \param m A VariantMatrix
     * \param refstate The ancestral state.
     * \return Fay and Wu's H, or nan if \a m is empty/invariant.
     *
     * \note It is an error to consider \a refstate as anything
     * other than the ancestral state for each site.
     *
     * This function is included via Sequence/summstats.hpp,
     * Sequence/summstats/classics.hpp
     * Sequence/summstats/thetah.hpp
     *
     * See \cite Fay2000-ef for details.
     * \ingroup popgenanalysis
     */
    double faywuh(const VariantMatrix& m, const std::int8_t refstate);

    /*! \brief Fay and Wu's H.
     * \param m A VariantMatrix
     * \param refstates The ancestral state at each site.
     * \return Fay and Wu's H, or nan if \a m is empty/invariant.
     *
     * \note It is an error to consider \a refstates as anything
     * other than the ancestral state at each site.
     *
     * This function is included via Sequence/summstats.hpp,
     * or Sequence/summstats/classics.hpp
     *
     * See \cite Fay2000-ef for details.
     * 
     * \ingroup popgenanalysis
     */
    double faywuh(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);

    /*!  \brief Calculate number of differences between all samples.
     * \param m A VariantMatrix
     * \return std::vector<std::int32_t>
     *
     * For \f$n\f$ samples in \a m, the output contains \f$n \choose 2\f$ elements.
     * More concretely, the elements are populated according to:
     *
     * \code
     * for(int i=0; i < m.nsam - 1; ++i)
     * {
     *      for(int j=i+1 ; j < m.nsam ; ++j)
     *      {
     *          //diffs[i] = number of diffs between sample i and j      
     *      }
     * }    
     * \endcode
     *  
     * Missing data to not contribute to differences between sequences.
     * Thus, low-quality data may lead to uninformative return values.
     *
     * Included via Sequence/summstats.hpp or 
     * Sequence/summstats/classics.hpp
     *
     * \ingroup popgenanalysis
     */
    std::vector<std::int32_t> difference_matrix(const VariantMatrix& m);

    /*! Returns whether or not haplotypes differ.
     * \param m A VariantMatrix
     * \return std::vector<std::int32_t>
     *
     * This function is conceptually indentical to difference_matrix,
     * but the return value contains 0 = identical, 1 = different 
     * rather than the actual number of differences. Thus, it is 
     * faster to calculate when the binary answer is needed.
     */
    std::vector<std::int32_t> is_different_matrix(const VariantMatrix& m);

    /*! \brief Assign a unique label to each haplotype.
     *
     * \param m A VariantMatrix
     * \return std::vector<std::int32_t>
     *
     * If there are \f$k\f$ unique samples in \a m, which represents
     * a sample of size \a nsam, the return value contains \a nsam
     * elements whose values are \f$[0,k)\f$ for "good" input.  Here,
     * "bad" input means that some of the samples consist entirely of 
     * missing data.  In that case, they are given the label of -1.
     *
     * This function is implemented via a call to difference_matrix.
     *
     * Included via Sequence/summstats.hpp or 
     * Sequence/summstats/classics.hpp
     *
     * \ingroup popgenanalysis
     */
    std::vector<std::int32_t> label_haplotypes(const VariantMatrix& m);

    /*! \brief Calculate the number of haplotypes in a sample.
     * \param m A VariantMatrix
     * \return The number of unique elements in the sample, std::int32_t
     *
     * This returns the number of unique columns in \a m.
     *
     * \note The value -1 is returned if \a m.nsam == 0.  If
     * \a m contains data, but all sites are invariant, the 
     * function returns 1.  See testClassicSummstatsEmptyVariantMatrix.cc
     * for examples.
     *
     * Include via Sequence/summstats.hpp or Sequence/summstats/classics.hpp
     *
     * See \cite Depaulis1998-ol for details.
     *
     * \ingroup popgenanalysis
     */
    std::int32_t number_of_haplotypes(const VariantMatrix& m);

    /*! \brief Calculate the haplotype diversity of a sample.
     * \param m A VariantMatrix
     * \return Haplotype heterozygosity, a double.
     *
     * The "haplotype heterozygosity" is calculated by counting
     * haplotype labels (see label_haplotypes).
     *
     * \note The value nan is returned if \a m.nsam == 0.  If
     * \a m contains data, but all sites are invariant, the 
     * function returns 0.0.  See testClassicSummstatsEmptyVariantMatrix.cc
     * for examples.
     *
     * Included via Sequence/summstats.hpp or 
     * Sequence/summstats/classics.hpp
     *
     * See \cite Depaulis1998-ol for details.
     *
     * \ingroup popgenanalysis
     */
    double haplotype_diversity(const VariantMatrix& m);

    /*! Hudson and Kaplan's Rmin statistic
     * \param m A VariantMatrix
     * \return Rmin, std::int32_t
     *
     * Included via Sequence/summstats.hpp or 
     * Sequence/summstats/classics.hpp
     *
     * See \cite Hudson1985-cq for details.
     * 
     * \ingroup popgenanalysis
     */
    std::int32_t rmin(const VariantMatrix& m);
} // namespace Sequence

#endif
