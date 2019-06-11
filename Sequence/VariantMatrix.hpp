#ifndef SEQUENCE_VARIANT_MATRIX_HPP__
#define SEQUENCE_VARIANT_MATRIX_HPP__

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <utility>
#include <vector>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include "VariantMatrixCapsule.hpp"

static_assert(sizeof(std::int8_t) == sizeof(char),
              "sizeof(char) is not 8 bits");

namespace Sequence
{
    /// \defgroup variantmatrix Variant Matrix
    /// \brief Representation and manipulation of variation data.
    /// \ingroup popgen

    class VariantMatrix
    /// \brief Matrix representation of variation data.
    ///
    /// The data structure is a row-major matrix.
    /// Variants are represented by 8-bit integers.
    /// Negative values represent missing data,
    /// thus allowing up to 128 non-missing states
    /// per variable site.
    ///
    /// Rows are sites, columns are either haplotypes
    /// or are populated as necessary in the case of
    /// genotype (unphased) data.  Storing variable
    /// sites in rows places the performance priority
    /// on sites over haplotypes.
    ///
    /// The only reserved character state is
    /// std::numeric_limits<std::int8_t>::min(),
    /// which is used to represent masked data.
    /// That reserved value is VariantMatrix::mask,
    /// which is a static constant.
    ///
    /// \ingroup variantmatrix
    /// \version 1.9.2
    {
      private:
        std::unique_ptr<VariantMatrixCapsule> capsule;
        std::int8_t
        set_max_allele(const std::int8_t max_allele_value)
        {
            if (max_allele_value < 0 && !capsule->empty())
                {
                    auto itr
                        = std::max_element(capsule->begin(), capsule->end());
                    return *itr;
                }
            //special case to allow construction of empty data sets
            //without throwing an exception
            else if (capsule->empty())

                {
                    return 0;
                }
            return max_allele_value;
        }
        std::int8_t max_allele_;

      public:
        /// Data stored in matrix form with rows as sites.
        //std::vector<std::int8_t> data;
        /// Position of sites.
        std::vector<double> positions;
        /// Number of sites in data set.
        std::size_t nsites;
        /// Sample size of data set.
        std::size_t nsam;
        /// Reserved value for masked data
        static const std::int8_t mask;
        /// The value type of the data.
        /// Helpful for generic programming
        using value_type = std::int8_t;
        template <typename data_input, typename positions_input>
        VariantMatrix(data_input&& data_, positions_input&& positions_,
                      const std::int8_t max_allele_value = -1)
            /// \brief "Perfect-forwarding" constructor.
            ///
            /// std::invalid_argument will be thrown if
            /// data.size() % positions.size() != 0.0.
            : capsule(new VectorCapsule(std::forward<data_input>(data_))),
              positions(std::forward<positions_input>(positions_)),
              nsites(positions.size()),
              nsam((nsites > 0) ? capsule->size() / positions.size() : 0),
              max_allele_{ set_max_allele(max_allele_value) }
        {
            if (max_allele() < 0)
                {
                    throw std::invalid_argument("max allele must be >= 0");
                }
            if ((!capsule->empty() && !positions.empty())
                && capsule->size() % positions.size() != 0.0)
                {
                    throw std::invalid_argument("incorrect dimensions");
                }
        }

        // Non range-checked access

        /// \brief Get data from marker `site` and haplotype `haplotype`.
        /// No range-checking is done.
        std::int8_t& get(const std::size_t site, const std::size_t haplotype);
        /// \brief Get data from marker `site` and haplotype `haplotype`.
        /// No range-checking is done.
        const std::int8_t& get(const std::size_t site,
                               const std::size_t haplotype) const;

        // Ranged-checked access after std::vector<T>::at.
        /// \brief Get data from marker `site` and haplotype `haplotype`.
        /// std::out_of_range is thrown if indexes are invalid.
        std::int8_t& at(const std::size_t site, const std::size_t haplotype);
        /// \brief Get data from marker `site` and haplotype `haplotype`.
        /// std::out_of_range is thrown if indexes are invalid.
        const std::int8_t& at(const std::size_t site,
                              const std::size_t haplotype) const;
        std::int8_t* data();
        const std::int8_t* data() const;
        bool empty() const;
        /// Max allelic value stored in matrix
        std::int8_t max_allele() const;

        void swap(VariantMatrix& rhs);
    };

    void swap(VariantMatrix& a, VariantMatrix& b);
} // namespace Sequence

#endif
