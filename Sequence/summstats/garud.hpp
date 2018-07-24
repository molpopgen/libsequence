#ifndef SEQUENCE_SUMMSTATS_GARUD_HPP
#define SEQUENCE_SUMMSTATS_GARUD_HPP

namespace Sequence
{
  struct GarudStats
  /*! 
    From http://arxiv.org/abs/1303.0906
    \note H1 = 1 - haplotype homozygosity, e.g. Depaulis and Veuille's "H"
    \ingroup popgenanalysis
  */
  {
    double H1,H12,H2H1;
    GarudStats();
    GarudStats(const double, const double, const double);
  };
}

#endif
