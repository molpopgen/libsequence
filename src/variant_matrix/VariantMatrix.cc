#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/StateCounts.hpp>
#include <stdexcept>

namespace Sequence
{
    std::vector<std::int32_t>
    AlleleCountMatrix::init_counts(const VariantMatrix& m)
    {
        if (m.max_allele < 0)
            {
                throw std::invalid_argument("matrix max_allele must be >= 0");
            }
        std::vector<std::int32_t> counts;
        counts.reserve(m.nsam * static_cast<std::size_t>(m.max_allele + 1));
        StateCounts c;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                c(r);
                for (std::size_t j = 0;
                     j < static_cast<std::size_t>(m.max_allele + 1); ++j)
                    {
                        counts.push_back(c.counts[j]);
                    }
            }
        return counts;
    }

    AlleleCountMatrix::AlleleCountMatrix(const VariantMatrix& m)
        : counts(init_counts(m)), row_size{
              static_cast<std::size_t>(m.max_allele) + 1
          }
    {
    }

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
