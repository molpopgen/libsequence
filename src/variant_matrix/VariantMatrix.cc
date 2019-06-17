#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/StateCounts.hpp>
#include <stdexcept>
#include <algorithm>

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

    const double&
    VariantMatrix::cposition(std::size_t i) const
    {
        auto cp = extract_const_ptr(this->pcapsule);
        return cp->operator[](i);
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

    const double*
    VariantMatrix::cpbegin() const
    {
        const auto cp = extract_const_ptr(pcapsule);
        return cp->begin();
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

    const double*
    VariantMatrix::cpend() const
    {
        const auto cp = extract_const_ptr(pcapsule);
        static_assert(std::is_const<decltype(cp)>::value, "foo");
        return cp->end();
    }

    // Non range-checked access
    std::int8_t&
    VariantMatrix::get(const std::size_t site, const std::size_t haplotype)
    {
        return capsule->operator()(site, haplotype);
    }

    const std::int8_t&
    VariantMatrix::get(const std::size_t site,
                       const std::size_t haplotype) const
    {
        return capsule->operator()(site, haplotype);
    }

    const std::int8_t&
    VariantMatrix::cget(const std::size_t site,
                        const std::size_t haplotype) const
    {
        auto cp = extract_const_ptr(this->capsule);
        return cp->operator()(site, haplotype);
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
        return capsule->operator()(site, haplotype);
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
        return capsule->operator()(site, haplotype);
    }

    const std::int8_t&
    VariantMatrix::cat(const std::size_t site,
                       const std::size_t haplotype) const
    {
        if (site >= nsites() || haplotype >= nsam())
            {
                throw std::out_of_range(
                    "VariantMatrix::at -- index out of range");
            }
        auto cp = extract_const_ptr(this->capsule);
        return cp->operator()(site, haplotype);
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

    VariantMatrix
    VariantMatrix::deepcopy() const
    {
        return VariantMatrix(capsule->clone(), pcapsule->clone(),
                             max_allele());
    }

    bool
    VariantMatrix::resizable() const
    {
        return capsule->resizable() && pcapsule->resizable();
    }

    void
    VariantMatrix::resize_capsules(bool remove_sites)
    {
        capsule->resize(remove_sites);
        pcapsule->resize(remove_sites);
    }

    std::size_t
    VariantMatrix::genotype_row_offset() const
    {
        return capsule->row_offset();
    }

    std::size_t
    VariantMatrix::genotype_col_offset() const
    {
        return capsule->col_offset();
    }

    std::size_t
    VariantMatrix::genotype_stride() const
    {
        return capsule->stride();
    }

    void
    swap(VariantMatrix& a, VariantMatrix& b)
    {
        a.swap(b);
    }

    bool
    operator==(const VariantMatrix& a, const VariantMatrix& b)
    {
        if (a.nsites() != b.nsites() || a.nsam() != b.nsam())
            return false;

        if (a.max_allele() != b.max_allele())
            return false;

        auto pdiff = std::mismatch(a.pbegin(), a.pend(), b.pbegin());
        if (pdiff.first != a.pend())
            return false;
        auto ddiff = std::mismatch(a.data(), a.data() + a.nsites() * a.nsam(),
                                   b.data());
        return ddiff.first == a.data() + a.nsites() * a.nsam();
    }

    bool
    operator!=(const VariantMatrix& a, const VariantMatrix& b)
    {
        return !(a == b);
    }

} // namespace Sequence
