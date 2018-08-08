#ifndef SEQUENCE_ALLELE_COUNT_MATRIX_HPP
#define SEQUENCE_ALLELE_COUNT_MATRIX_HPP

#include <cstdint>
#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    class AlleleCountMatrix
    /// \brief Matrix representation of allele counts in a VariantMatrix
    /// To be constructed
    {
      private:
        static std::vector<std::int32_t> init_counts(const VariantMatrix& m);

      public:
        const std::vector<std::int32_t> counts;
        const std::size_t ncol;
        const std::size_t nrow;
        const std::size_t nsam;
        explicit AlleleCountMatrix(const VariantMatrix& m);
    };
}

#endif
