#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/StateCounts.hpp>

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
        : counts(init_counts(m)),
          row_size(static_cast<std::size_t>(m.max_allele) + 1), nsam(m.nsam)
    {
    }
}
