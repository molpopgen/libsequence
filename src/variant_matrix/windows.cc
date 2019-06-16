#include <Sequence/NonOwningCapsules.hpp>
#include <Sequence/variant_matrix/windows.hpp>

namespace Sequence
{
    VariantMatrix
    make_window(const VariantMatrix& m, const double beg, const double end)
    {
        return make_slice(m, beg, end, 0, m.nsam());
    }

    VariantMatrix
    make_slice(const VariantMatrix& m, const double beg, const double end,
               const std::size_t i, const std::size_t j)
    {
        if (end < beg)
            {
                throw std::invalid_argument("end must be >= beg");
            }
        if (!(j > i))
            {
                throw std::invalid_argument("i must be < j");
            }
        if (j > m.nsam())
            {
                throw std::invalid_argument("slice indexes out of range");
            }
        auto pb = std::lower_bound(m.pbegin(), m.pend(), beg);
        auto pe = std::upper_bound(pb, m.pend(), end);
        if (pb == m.pend())
            {
                std::unique_ptr<GenotypeCapsule> gc(
                    new NonOwningGenotypeCapsule(m.cdata(), 0, 0, 0, 0, 0));
                std::unique_ptr<PositionCapsule> pc(
                    new NonOwningPositionCapsule(pb, 0));
                return VariantMatrix(std::move(gc), std::move(pc), -1);
            }
        std::size_t nsites = pe - pb;
        std::size_t nsam = j - i;
        std::size_t row_offset = pb - m.pbegin();
        std::unique_ptr<GenotypeCapsule> gc(new NonOwningGenotypeCapsule(
            m.cdata(), nsites, nsam, row_offset, i, m.nsam()));
        std::unique_ptr<PositionCapsule> pc(
            new NonOwningPositionCapsule(pb, pe - pb));
        return VariantMatrix(std::move(gc), std::move(pc), m.max_allele());
    }
} // namespace Sequence
