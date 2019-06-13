#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/StateCounts.hpp>
#include <stdexcept>

namespace Sequence
{
    const std::int8_t VariantMatrix::mask
        = std::numeric_limits<std::int8_t>::min();

    std::int8_t
    VariantMatrix::max_allele() const
    {
        return max_allele_;
    }

    std::size_t
    VariantMatrix::nsites() const
    {
        return capsule->nsites();
    }
    std::size_t&
    VariantMatrix::nsites()
    {
        return capsule->nsites();
    }

    std::size_t
    VariantMatrix::nsam() const
    {
        return capsule->nsam();
    }

    std::size_t&
    VariantMatrix::nsam()
    {
        return capsule->nsam();
    }

    double
    VariantMatrix::position(std::size_t i) const
    {
        return pcapsule->operator[](i);
    }

    double&
    VariantMatrix::position(std::size_t i)
    {
        return pcapsule->operator[](i);
    }

    double*
    VariantMatrix::pbegin()
    {
        return pcapsule->begin();
    }

    const double*
    VariantMatrix::pbegin() const
    {
        return pcapsule->begin();
    }

    double*
    VariantMatrix::pend()
    {
        return pcapsule->end();
    }

    const double*
    VariantMatrix::pend() const
    {
        return pcapsule->end();
    }

    // Non range-checked access
    std::int8_t&
    VariantMatrix::get(const std::size_t site, const std::size_t haplotype)
    {
        return capsule->operator[](site* nsam() + haplotype);
    }

    const std::int8_t&
    VariantMatrix::get(const std::size_t site,
                       const std::size_t haplotype) const
    {
        return capsule->operator[](site* nsam() + haplotype);
    }

    // Ranged-checked access after std::vector<T>::at.
    std::int8_t&
    VariantMatrix::at(const std::size_t site, const std::size_t haplotype)
    {
        if (site >= nsites() || haplotype >= nsam())
            {
                throw std::out_of_range(
                    "VariantMatrix::at -- index out of range");
            }
        return capsule->operator[](site* nsam() + haplotype);
    }

    const std::int8_t&
    VariantMatrix::at(const std::size_t site,
                      const std::size_t haplotype) const
    {
        if (site >= nsites() || haplotype >= nsam())
            {
                throw std::out_of_range(
                    "VariantMatrix::at -- index out of range");
            }
        return capsule->operator[](site* nsam() + haplotype);
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

    const std::int8_t*
    VariantMatrix::cdata() const
    {
        return capsule->cdata();
    }

    bool
    VariantMatrix::empty() const
    {
        return capsule->empty();
    }

    void
    VariantMatrix::swap(VariantMatrix& rhs)
    {
        this->capsule.swap(rhs.capsule);
        this->pcapsule.swap(rhs.pcapsule);
        std::swap(this->max_allele_, rhs.max_allele_);
    }

    bool
    VariantMatrix::resizable() const
    {
        return capsule->resizable() && pcapsule->resizable();
    }

    void
    VariantMatrix::resize_capsules()
    {
        capsule->resize();
        pcapsule->resize();
    }

    void
    swap(VariantMatrix& a, VariantMatrix& b)
    {
        a.swap(b);
    }
} // namespace Sequence
