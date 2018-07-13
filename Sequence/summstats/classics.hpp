#ifndef SEQUENCE_SUMMSTATS_CLASSICS_HPP__
#define SEQUENCE_SUMMSTATS_CLASSICS_HPP__

#include <Sequence/VariantMatrix.hpp>

#include "thetapi.hpp"
#include "thetaw.hpp"
#include "thetah.hpp"
#include "nvariablesites.hpp"
#include "allele_counts.hpp"

namespace Sequence
{
    double tajd(const VariantMatrix&);

    double faywuh(const VariantMatrix& m, const std::int8_t refstate);

    double faywuh(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);
} // namespace Sequence

#endif
