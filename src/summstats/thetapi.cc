#include <unordered_map>
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
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_ConstRowView(m, site);
                StateCounts counts(site_view);
                //StateCounts will collapse all missing data
                //values into a single value -1.  Instead of
                //checking if each value is not missing, we
                //simply remove all missing data:
                counts.counts.erase(-1);
                double sshom = 0.0; //sum of site homozygosity
                for (auto& x : counts.counts)
                    {
                        sshom += static_cast<double>(x.second)
                                 * static_cast<double>(x.second - 1.0);
                    }
                pi += 1.0
                      - sshom
                            / (static_cast<double>(counts.n)
                               * (static_cast<double>(counts.n) - 1.0));
            }
        return pi;
    }

} // namespace Sequence
