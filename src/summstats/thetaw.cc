#include <algorithm>
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
        StateCounts counts;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                counts(site_view);
                if (counts.counts.size() > 1)
                    {
                        auto nstates = counts.counts.size()
                                       - std::count(counts.counts.begin(),
                                                    counts.counts.end(), 0);
                        auto denom = summstats_aux::a_sub_n(counts.n);
                        w += static_cast<double>(nstates - 1) / denom;
                    }
            }
        return w;
    }

    double
    thetaw(const AlleleCountMatrix& ac)
    {
        double w = 0.0;
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.row_size)
            {
                std::int32_t nsam = 0, nstates = 0;
                for (std::size_t j = i; j < i + ac.row_size; ++j)
                    {
                        if (ac.counts[j] > 0)
                            {
                                nsam += ac.counts[j];
                                nstates++;
                            }
                    }
                if (nstates > 1)
                    {
                        auto denom = summstats_aux::a_sub_n(nsam);
                        w += static_cast<double>(nstates - 1) / denom;
                    }
            }
        return w;
    }

} // namespace Sequence
