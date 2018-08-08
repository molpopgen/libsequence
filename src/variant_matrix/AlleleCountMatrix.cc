#include <stdexcept>
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
                if (static_cast<std::int8_t>(c.max_allele_idx) > m.max_allele)
                    {
                        throw std::runtime_error("found allele value greater "
                                                 "than matrix.max_allele");
                    }
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
          ncol(!m.data.empty() ? static_cast<std::size_t>(m.max_allele) + 1
                               : 0),
          nrow(!m.data.empty() ? m.data.size() / ncol : 0), nsam(m.nsam)
    {
    }

    std::pair<std::vector<std::int32_t>::const_iterator,
              std::vector<std::int32_t>::const_iterator>
    AlleleCountMatrix::row(const std::size_t i) const
    {
        if (i >= nrow)
            {
                throw std::out_of_range("row index out of range");
            }
        return std::make_pair(counts.begin() + i * ncol,
                              counts.begin() + i * ncol + ncol);
    }
} // namespace Sequence
