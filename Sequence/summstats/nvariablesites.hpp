#ifndef SEQUENCE_SUMMSTATS_NVARIABLESITES_HPP__
#define SEQUENCE_SUMMSTATS_NVARIABLESITES_HPP__

#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    std::uint32_t nvariable_sites(const VariantMatrix& m);
    std::uint32_t nbiallelic_sites(const VariantMatrix& m);
    std::uint32_t total_number_of_mutations(const VariantMatrix& m);
} // namespace Sequence

#endif
