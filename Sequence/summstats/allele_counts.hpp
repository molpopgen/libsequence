#ifndef SEQUENCE_SUMMSTATS_ALLELE_COUNTS_HPP__
#define SEQUENCE_SUMMSTATS_ALLELE_COUNTS_HPP__

#include <vector>
#include <utility>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
	struct AlleleCounts
	{
		int nstates, nmissing;
	};

    std::vector<AlleleCounts>
    allele_counts(const VariantMatrix& m);

    std::vector<AlleleCounts>
    non_reference_allele_counts(const VariantMatrix& m, const std::int8_t refstate);

    std::vector<AlleleCounts>
    non_reference_allele_counts(const VariantMatrix& m, const std::vector<std::int8_t> & refstates);
} // namespace Sequence
#endif
