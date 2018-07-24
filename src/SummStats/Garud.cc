#include <algorithm>
#include <set>
#include <string>
#include <cmath>
#include <numeric>
#include <Sequence/SummStats/Garud.hpp>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

using namespace std;

namespace Sequence
{
    GarudStats
    H1H12(const SimData &d)
    /*!
      H1 is total haplotype homozygosity.
      H2 is haplotype homozygosity, combining two most common haplotypes.  H2 =
      H1 + 2p1p2
      H2H1 = H2/H1, where H2 is haplotype homozygosity for all but most common
      haplotype.
      H2H1 = (H1 - p1^2)/H1
     */
    {
        if (d.empty())
            return GarudStats();
        set<string> uhaps(d.begin(), d.end());
        vector<string> vuhaps(uhaps.size());
        vector<double> hapcounts(uhaps.size());
        std::move(uhaps.begin(), uhaps.end(), vuhaps.begin());
        tbb::parallel_for(
            tbb::blocked_range<std::size_t>(0, uhaps.size()),
            [&d, &vuhaps,
             &hapcounts](const tbb::blocked_range<std::size_t> &r) {
                for (auto i = r.begin(); i < r.end(); ++i)
                    {
                        hapcounts[i] = static_cast<double>(
                            std::count(d.begin(), d.end(), vuhaps[i]));
                    }
            });
        const double denom = static_cast<double>(d.size() * (d.size() - 1));
        double H1 = tbb::parallel_reduce(
            tbb::blocked_range<double *>(hapcounts.data(),
                                         hapcounts.data() + hapcounts.size()),
            0.,
            [denom](const tbb::blocked_range<double *> &r,
                    double value) -> double {
                return std::accumulate(r.begin(), r.end(), value,
                                       [denom](double s, double c) {
                                           return s + c * (c - 1.) / denom;
                                       });
            },
            std::plus<double>());
        // GarudStats G;
        sort(hapcounts.begin(), hapcounts.end(),
             std::bind(greater<double>(), std::placeholders::_1,
                       std::placeholders::_2));
        double H12 = H1
                     + 2. * hapcounts[0] * hapcounts[1]
                           / std::pow(double(d.size()), 2.);
        double H2H1 = (H1
                       - double(hapcounts[0] * (hapcounts[0] - 1))
                             / double(d.size() * (d.size() - 1)))
                      / H1;
        return GarudStats(H1, H12, H2H1);
    }
}
