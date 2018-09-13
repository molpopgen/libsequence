#ifndef SEQUENCE_ALLELE_COUNT_MATRIX_HPP
#define SEQUENCE_ALLELE_COUNT_MATRIX_HPP

#include <cstdint>
#include <vector>
#include <utility>
#include <stdexcept>
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
        using value_type = std::vector<std::int32_t>::value_type;
        const std::size_t ncol;
        const std::size_t nrow;
        const std::size_t nsam;
        explicit AlleleCountMatrix(const VariantMatrix& m);

        /// This constructor is for advanced use only,
        /// such as constructing from a slice of a
        /// pre-existing AlleleCountMatrix.
        template <typename T>
        AlleleCountMatrix(T&& t, const std::size_t nc_, const std::size_t nr_,
                          const std::size_t n_)
            : counts(std::forward<T>(t)), ncol{ nc_ }, nrow{ nr_ }, nsam{ n_ }
        {
            if (ncol * nrow != counts.size())
                {
                    throw std::invalid_argument(
                        "incorrect dimensions for AlleleCountMatrix");
                }
        }
        std::pair<std::vector<std::int32_t>::const_iterator,
                  std::vector<std::int32_t>::const_iterator>
        row(const std::size_t) const;
    };
} // namespace Sequence

#endif
