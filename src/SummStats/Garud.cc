#include <algorithm>
#include <set>
#include <string>
#include <cmath>
#include <numeric>
#include <Sequence/SummStats/Garud.hpp>

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
        std::move(uhaps.begin(), uhaps.end(), vuhaps.begin());
        vector<double> hapcounts;
		hapcounts.reserve(uhaps.size());
		for(auto & uh : uhaps)
		{
			hapcounts.push_back(static_cast<double>(std::count(d.begin(),d.end(),uh)));
		}
        const double denom = static_cast<double>(d.size() * (d.size() - 1));
        double H1 = 0.0;
		for(auto c : hapcounts)
		{
			H1 += c*(c-1.0);
		}
		H1 /= denom;
        
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
