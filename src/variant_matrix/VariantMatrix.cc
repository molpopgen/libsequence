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
        return capsule->operator[](site* nsam + haplotype);
    }

    const std::int8_t&
    VariantMatrix::get(const std::size_t site,
                       const std::size_t haplotype) const
    {
        return capsule->operator[](site* nsam + haplotype);
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
        return capsule->operator[](site* nsam + haplotype);
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
        return capsule->operator[](site* nsam + haplotype);
    }

    std::int8_t*
    VariantMatrix::data()
    {
        return capsule->data();
    }

    const std::int8_t*
    VariantMatrix::data() const
    {
        return capsule->data();
    }

    bool
    VariantMatrix::empty() const
    {
        return capsule->empty();
    }
} // namespace Sequence
