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
                for (auto& x : counts.counts)
                    {
                        sshom += static_cast<double>(x)
                                 * static_cast<double>(x - 1.0);
                    }
                double nnm1 = (static_cast<double>(counts.n)
                               * (static_cast<double>(counts.n) - 1.0));
                pi += 1.0 - sshom / nnm1;
            }
        return pi;
    }

} // namespace Sequence
