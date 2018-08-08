#include <cmath>
#include <limits>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    tajd(const AlleleCountMatrix& ac)
    {
        double pi = 0.0;
        int S = 0;
        std::int32_t max_nsam = 0;
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.ncol)
            {
                std::int32_t nsam = 0;
                double homozygosity = 0.0;
                int nstates = 0;
                for (std::size_t j = i; j < i + ac.ncol; ++j)
                    {
                        if (ac.counts[j] > 0)
                            {
                                ++nstates;
                                nsam += ac.counts[j];
                                homozygosity += static_cast<double>(
                                    ac.counts[j] * (ac.counts[j] - 1));
                            }
                    }

                if (nstates)
                    {
                        max_nsam = std::max(max_nsam, nsam);
                        S += nstates - 1;
                        pi += 1.0
                              - homozygosity
                                    / static_cast<double>(nsam * (nsam - 1));
                    }
            }
        if (!S)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto a1 = summstats_aux::a_sub_n(static_cast<std::uint32_t>(max_nsam));
        double w = static_cast<double>(S) / a1;
        auto a2 = summstats_aux::b_sub_n(static_cast<std::uint32_t>(max_nsam));
        auto dn = static_cast<double>(max_nsam);
        double b1 = (dn + 1.0) / (3.0 * (dn - 1.0));
        double b2
            = (2.0 * (std::pow(dn, 2.0) + dn + 3.0)) / (9.0 * dn * (dn - 1.0));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (dn + 2.0) / (a1 * dn) + a2 / std::pow(a1, 2.0);
        double e1 = c1 / a1;
        double e2 = c2 / (std::pow(a1, 2.0) + a2);
        double denominator = std::pow((e1 * S + e2 * S * (S - 1.0)), 0.5);
        return (pi - w) / denominator;
    }
} // namespace Sequence
