#ifndef __SEQUENCE_GARUD_HPP__
#define __SEQUENCE_GARUD_HPP__

#include <Sequence/SimData.hpp>

namespace Sequence
{
  struct GarudStats
  /*! 
    From http://arxiv.org/abs/1303.0906
    \note H1 = 1 - haplotype homozygosity, e.g. Depaulis and Veuille's "H"
  */
  {
    mutable double H1,H12,H2H1;
    GarudStats();
    GarudStats(const double, const double, const double);
  };

  /*
    Garud et al. DOI: 10.1371/journal.pgen.1005004
    Messer & Petrov DOI: 10.1016/j.tree.2013.08.003
    Note that H1 = 1 - haplotype homozygosity, e.g. Depaulis and Veuille's "H"

    \return An object of type Sequence::GarudStats
  */
  GarudStats H1H12(const SimData & d);
}

#endif
