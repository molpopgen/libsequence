/// \file Sequence/summstats/garud.hpp
/// \brief H1, H12, and H2/H1 stats
#ifndef SEQUENCE_SUMMSTATS_GARUD_HPP
#define SEQUENCE_SUMMSTATS_GARUD_HPP

#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
  struct GarudStats
  /*! 
    Statistics from \cite Garud2015-ob
    \note H1 = 1 - haplotype homozygosity, e.g. "H" from \cite Depaulis1998-ol
    \ingroup popgenanalysis
  */
  {
    double H1,H12,H2H1;
    GarudStats();
    GarudStats(const double, const double, const double);
  };

  /*! \brief Calculate H1, H12, and H2/H1
   * \param m A VariantMatrix
   * \return GarudStats
   *
   * See \cite Garud2015-ob for details.
   */
  GarudStats garud_statistics(const VariantMatrix & m);
}

#endif
