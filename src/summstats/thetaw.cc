#include <stdexcept>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    thetaw(const VariantMatrix& m)
    {
        double w = 0.0;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                StateCounts counts(site_view);
                counts.counts.erase(-1);
                if (counts.counts.size() > 1)
                    {
                        auto nstates = counts.counts.size();
                        auto denom = summstats_aux::a_sub_n(counts.n);
                        w += static_cast<double>(nstates - 1) / denom;
                    }
            }
        return w;
    }

} // namespace Sequence
