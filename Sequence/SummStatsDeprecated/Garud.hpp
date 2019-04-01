#ifndef __SEQUENCE_GARUD_HPP__
#define __SEQUENCE_GARUD_HPP__

#include <Sequence/SimData.hpp>
#include <Sequence/summstats/garud.hpp>

namespace Sequence
{
  /*
    Garud et al. DOI: 10.1371/journal.pgen.1005004
    Messer & Petrov DOI: 10.1016/j.tree.2013.08.003
    Note that H1 = 1 - haplotype homozygosity, e.g. Depaulis and Veuille's "H"
    \ingroup popgenanalysis
    \return An object of type Sequence::GarudStats
  */
  GarudStats H1H12(const SimData & d)__attribute__ ((deprecated));
}

#endif
