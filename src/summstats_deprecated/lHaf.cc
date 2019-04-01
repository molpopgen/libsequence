#include <Sequence/SummStatsDeprecated/lHaf.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Sequence
{
    std::vector<double>
    lHaf(const SimData &data, const double l)
    {
        //using range_type = tbb::blocked_range<SimData::const_site_iterator>;
        //using data_range_type
        //    = tbb::blocked_range<std::vector<std::string>::const_iterator>;
        // Get derived mutation frequency counts per site
        std::vector<unsigned> dcounts;
        dcounts.reserve(data.numsites());
        for (auto i = data.sbegin(); i < data.send(); ++i)
            {
                dcounts.push_back(static_cast<unsigned>(
                    std::count(i->second.begin(), i->second.end(), '1')));
            }
        // Get the values for each element in the data
        std::vector<double> rv;
        rv.reserve(data.size());
        for (auto &i : data)
            {
                auto j
                    = std::find_if(i.cbegin(), i.cend(),
                                   [](const char &ch) { return ch == '1'; });
                double score = 0.0;
                while (j != i.cend())
                    {
                        size_t d2 = size_t(j - i.cbegin());
                        score += std::pow(static_cast<double>(dcounts[d2]), l);
                        j = std::find(j + 1, i.cend(), '1');
                    }
                rv.push_back(score);
            }
        return rv;
    }
} // namespace Sequence
