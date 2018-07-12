#ifndef SEQUENCE_VARIANT_MATRIX_HPP__
#define SEQUENCE_VARIANT_MATRIX_HPP__

#include <cstdint>
#include <cstddef>
#include <utility>
#include <vector>
#include <limits>
#include <stdexcept>
#include <type_traits>

static_assert(sizeof(std::int8_t) == sizeof(char),
              "sizeof(char) is not 8 bits");

namespace Sequence
{
    /// \defgroup variantmatrix Variant Matrix
    /// \brief Types and functions related to manipulation of variation data

    struct VariantMatrix
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
        /// Data stored in matrix form with rows as sites.
        std::vector<std::int8_t> data;
        /// Position of sites.
        std::vector<double> positions;
        /// Number of sites in data set.
        std::size_t nsites;
        /// Sample size of data set.
        std::size_t nsam;
        /// Reserved value for masked data
        static const std::int8_t mask;

        template <typename data_input, typename positions_input>
        VariantMatrix(data_input&& data_, positions_input&& positions_)
            /// \brief "Perfect-forwarding" constructor.
            ///
            /// std::invalid_argument will be thrown if
            /// data.size() % positions.size() != 0.0.
            : data(std::forward<data_input>(data_)),
              positions(std::forward<positions_input>(positions_)),
              nsites(positions.size()), nsam(data.size() / positions.size())
        {
            if (data.size() % positions.size() != 0.0)
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
    };
}

#endif
