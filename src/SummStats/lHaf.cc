#include <Sequence/SummStats/lHaf.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
/*
  TODO:  A thread_site and thread_hap wrapper that apply some fxn to a site
  or a hap, using <thread>, filling up a vector of return values.

  Need to read/remind myself of this:

  http://www.aristeia.com/TalkNotes/ACCU2011_CPUCaches.pdf

  This sounds cool, too: c++ concurrency book
 */
namespace Sequence
{
    std::vector<double>
    lHaf(const SimData &data, const double l)
    {
        using range_type = tbb::blocked_range<SimData::const_site_iterator>;
        using data_range_type
            = tbb::blocked_range<std::vector<std::string>::const_iterator>;
        // Get derived mutation frequency counts per site
        std::vector<unsigned> dcounts(data.numsites());
        tbb::parallel_for(
            range_type(data.sbegin(), data.send()),
            [&dcounts, &data](const range_type &r) {
                for (auto i = r.begin(); i < r.end(); ++i)
                    {
                        dcounts[static_cast<std::size_t>(
                            std::distance(data.sbegin(), i))]
                            = static_cast<unsigned>(std::count(
                                i->second.begin(), i->second.end(), '1'));
                    }
            });
        // Get the values for each element in the data
        std::vector<double> rv(data.size(), 0.);
        auto dbegin = data.begin();
        tbb::parallel_for(
            data_range_type(data.begin(), data.end()),
            [dbegin, l,&rv, &dcounts](const data_range_type &r) {
                for (auto i = r.begin(); i < r.end(); ++i)
                    {
                        auto d = static_cast<std::size_t>(
                            std::distance(dbegin, i));
                        auto j = std::find_if(
                            i->cbegin(), i->cend(),
                            [](const char &ch) { return ch == '1'; });
                        while (j != i->cend())
                            {
                                size_t d2 = size_t(j - i->cbegin());
                                rv[d] += std::pow(double(dcounts[d2]), l);
                                j = std::find_if(
                                    j + 1, i->cend(),
                                    [](const char &ch) { return ch == '1'; });
                            }
                    }
            });

        return rv;
    }
}
