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

    std::vector<std::int32_t> difference_matrix(const VariantMatrix& m);
    std::vector<std::int32_t> label_haplotypes(const VariantMatrix& m);
    std::uint32_t number_of_haplotypes(const VariantMatrix& m);
    double haplotype_diversity(const VariantMatrix& m);
} // namespace Sequence

#endif
