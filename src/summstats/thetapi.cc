#include <stdexcept>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    double
    thetapi(const VariantMatrix& m)
    {
        double pi = 0.0;
        StateCounts counts;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_ConstRowView(m, site);
                counts(site_view);
                double sshom = 0.0; //sum of site homozygosity
                for (std::size_t i = 0; i < counts.max_allele_idx + 1; ++i)
                    {
                        auto x = counts.counts[i];
                        sshom += static_cast<double>(x)
                                 * static_cast<double>(x - 1.0);
                    }
                double nnm1 = (static_cast<double>(counts.n)
                               * (static_cast<double>(counts.n) - 1.0));
                pi += 1.0 - sshom / nnm1;
            }
        return pi;
    }

    double
    thetapi(const AlleleCountMatrix& ac)
    {
        double pi = 0.0;
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.row_size)
            {
                std::int32_t nsam = 0;
                double homozygosity = 0.0;
                for (std::size_t j = i; j < i + ac.row_size; ++j)
                    {
                        nsam += ac.counts[j];
                        homozygosity += static_cast<double>(
                            ac.counts[j] * (ac.counts[j] - 1));
                    }

                pi += 1.0
                      - homozygosity / static_cast<double>(nsam * (nsam - 1));
            }
        return pi;
    }

} // namespace Sequence
