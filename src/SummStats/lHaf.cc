#include <Sequence/SummStats/lHaf.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>

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
        // Get derived mutation frequency counts per site
        std::vector<unsigned> dcounts;
        std::for_each(
            data.sbegin(), data.send(), [&dcounts](const polymorphicSite &p) {
                dcounts.push_back(unsigned(
                    std::count(p.second.begin(), p.second.end(), '1')));
            });
        // Get the values for each element in the data
        std::vector<double> rv(data.size(), 0.);
        for (unsigned i = 0; i < data.size(); ++i)
            {
                auto j
                    = std::find_if(data[i].cbegin(), data[i].cend(),
                                   [](const char &ch) { return ch == '1'; });
                while (j != data[i].cend())
                    {
                        size_t d = size_t(j - data[i].cbegin());
                        rv[i] += std::pow(double(dcounts[d]), l);
                        j = std::find_if(
                            j + 1, data[i].cend(),
                            [](const char &ch) { return ch == '1'; });
                    }
            }
        return rv;
    }
}
