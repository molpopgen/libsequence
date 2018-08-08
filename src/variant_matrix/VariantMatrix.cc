#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/StateCounts.hpp>
#include <stdexcept>

namespace Sequence
{
    const std::int8_t VariantMatrix::mask
        = std::numeric_limits<std::int8_t>::min();
    // Non range-checked access
    std::int8_t&
    VariantMatrix::get(const std::size_t site, const std::size_t haplotype)
    {
        return data[site * nsam + haplotype];
    }

    const std::int8_t&
    VariantMatrix::get(const std::size_t site,
                       const std::size_t haplotype) const
    {
        return data[site * nsam + haplotype];
    }

    // Ranged-checked access after std::vector<T>::at.
    std::int8_t&
    VariantMatrix::at(const std::size_t site, const std::size_t haplotype)
    {
        if (site >= nsites || haplotype >= nsam)
            {
                throw std::out_of_range(
                    "VariantMatrix::at -- index out of range");
            }
        return data.at(site * nsam + haplotype);
    }

    const std::int8_t&
    VariantMatrix::at(const std::size_t site,
                      const std::size_t haplotype) const
    {
        if (site >= nsites || haplotype >= nsam)
            {
                throw std::out_of_range(
                    "VariantMatrix::at -- index out of range");
            }
        return data.at(site * nsam + haplotype);
    }

    //AlleleCountMatrix
    //VariantMatrix::count_alleles() const
    //{
    //    return AlleleCountMatrix(std::move(counts),
    //                             static_cast<std::size_t>(max_allele + 1));
    //}
} // namespace Sequence
