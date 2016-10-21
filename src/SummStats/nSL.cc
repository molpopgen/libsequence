#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>
// For parallelizing nSL:
#include <functional>
#include <iostream>
#include <tbb/parallel_for.h>

using namespace std;

namespace
{
    inline double
    maxabs(double score, double mean, double sd, double rv)
    {
        if (isfinite(score))
            {
                double zscore = (score - mean) / sd;
                if (!isfinite(rv) || fabs(zscore) > fabs(rv))
                    return zscore;
            }
        return rv;
    }

    double
    update_s2(std::string::const_iterator start,
              std::string::const_iterator left,
              std::string::const_iterator right, const Sequence::SimData &d,
              const std::unordered_map<double, double> &gmap)
    {
        auto p1 = d.position(
            std::vector<double>::size_type(distance(start, right)) - 1);
        auto p2 = d.position(
            std::vector<double>::size_type(distance(start, left)));
        if (gmap.empty())
            {
                // return phyisical distance
                return fabs(p1 - p2);
            }
        // return distance along genetic map,
        // in whatever units those are.
        auto fp1 = gmap.find(p1);
        auto fp2 = gmap.find(p2);
        if (fp1 == gmap.end() || fp2 == gmap.end())
            {
                throw std::runtime_error(
                    "position could not be found in genetic map, "
                    + std::string(__FILE__) + " line "
                    + std::to_string(__LINE__));
            }
        return fabs(fp1->second - fp2->second);
    }
    /*
      Mechanics of the nSL statistic

      RV = nSL,iHS, as defined in doi:10.1093/molbev/msu077
    */
    std::array<double, 4>
    __nSLdetails(const std::size_t &core, const Sequence::SimData &d,
                 // const vector<size_t> &coretype,
                 const std::unordered_map<double, double> &gmap)
    {
        auto csize = d.size();
        // This tracks s,s2 for ancestral and derived
        // mutation, resp:
        std::array<double, 4> rv = { 0., 0., 0., 0. };
        // number of comparisons for ancestral and
        // derived, resp:
        std::array<unsigned, 2> nc = { 0u, 0u };
        for (size_t i = 0; i < csize; ++i)
            {
                auto start = d[i].cbegin();
                auto bi
                    = start
                      + static_cast<decltype(start)::difference_type>(core);
                size_t iIsDer = (*bi == '1');
                for (size_t j = i + 1; j < csize; ++j)
                    {
                        auto bj
                            = d[j].cbegin()
                              + static_cast<decltype(start)::difference_type>(
                                    core);
                        size_t jIsDer = (*bj == '1');
                        if (iIsDer == jIsDer)
                            {
                                auto eri = d[i].crend();
                                auto ei = d[i].cend();
                                auto right = mismatch(bi, ei, bj);
                                string::const_reverse_iterator ri1(bi),
                                    ri2(bj);
                                auto left = mismatch(ri1, eri, ri2);
                                if (left.first != eri && right.first != ei)
                                    {
                                        rv[2 * iIsDer] += static_cast<double>(
                                            distance(left.first.base(),
                                                     right.first)
                                            + 1);
                                        rv[2 * iIsDer + 1] += update_s2(
                                            start, left.first.base(),
                                            right.first, d, gmap);
                                        nc[iIsDer]++;
                                    }
                            }
                    }
            }
        rv[0] /= static_cast<double>(nc[0]);
        rv[1] /= static_cast<double>(nc[0]);
        rv[2] /= static_cast<double>(nc[1]);
        rv[3] /= static_cast<double>(nc[1]);
        return rv;
    }
}

namespace Sequence
{
    /*
      The nSL statistic of doi:10.1093/molbev/msu077
    */
    pair<double, double>
    nSL(const std::size_t &core, const SimData &d,
        const std::unordered_map<double, double> &gmap)
    {
        auto nsl = __nSLdetails(core, d, gmap);
        return make_pair(log(nsl[0]) - log(nsl[2]), log(nsl[1]) - log(nsl[3]));
    }

    vector<pair<double, double>>
    nSL_t(const SimData &d, const std::unordered_map<double, double> &gmap)
    {
        vector<pair<double, double>> rv(d.numsites());
        tbb::parallel_for(
            tbb::blocked_range<std::size_t>(0, rv.size()),
            [&rv, &d, &gmap](const tbb::blocked_range<std::size_t> &r) {
                for (std::size_t i = r.begin(); i < r.end(); ++i)
                    {
                        auto temp = __nSLdetails(i, d, gmap);
                        rv[i] = make_pair(log(temp[0]) - log(temp[2]),
                                          log(temp[1]) - log(temp[3]));
                    }
            });
        return rv;
    }
    /*
      Return max. abs value of standardized nSL and iHS, with the latter as
      defined by Ferrer-Admetella et al.
    */
    pair<double, double>
    snSL(const SimData &d, const double minfreq, const double binsize,
         const std::unordered_map<double, double> &gmap)
    {
        if (d.empty())
            return make_pair(std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN());
        vector<polymorphicSite> filtered;
        vector<unsigned> dcounts;
        dcounts.reserve(d.numsites());
        for_each(d.sbegin(), d.send(), [&](const polymorphicSite &p) {
            unsigned dcount = static_cast<unsigned>(
                count(p.second.begin(), p.second.end(), '1'));
            if (dcount && dcount < d.size())
                {
                    double f = static_cast<double>(dcount)
                               / static_cast<double>(d.size());
                    if (min(f, 1. - f) >= minfreq)
                        {
                            filtered.push_back(p);
                            dcounts.push_back(dcount);
                        }
                }
        });
        if (filtered.empty())
            return make_pair(std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN());
        SimData __filtered(filtered.begin(), filtered.end());
        // Get the stats
        auto nSLstats = nSL_t(__filtered, gmap);
        // Associate the stats with their DAFs
        using pp = pair<double, pair<double, double>>;
        vector<pp> binning;
        for (unsigned i = 0; i < __filtered.numsites(); ++i)
            {
                // pair<double, double> rvi = nSL(i, __filtered, gmap);
                binning.push_back(
                    make_pair(static_cast<double>(dcounts[i])
                                  / static_cast<double>(__filtered.size()),
                              nSLstats[i]));
            }
        // sort based on DAF
        std::sort(binning.begin(), binning.end(),
                  [](const pp &a, const pp &b) { return a.first < b.first; });
        double rv = std::numeric_limits<double>::quiet_NaN(),
               rv2 = std::numeric_limits<double>::quiet_NaN();
        // Now, bin, standardise, and move on...
        vector<pp> thisbin;
        auto bstart = binning.begin();
        size_t ttlSNPs = 0;
        for (double l = minfreq; l < 1.; l += binsize)
            {
                thisbin.clear();
                auto first = lower_bound(
                    bstart, binning.end(), l,
                    [](const pp &p, const double d) { return p.first < d; });
                auto last = upper_bound(
                    first, binning.end(), l + binsize,
                    [](const double d, const pp &p) { return d <= p.first; });
                thisbin.insert(thisbin.end(), first, last);
                ttlSNPs += thisbin.size();
                bstart = last;
                if (thisbin.size() > 1) // otherwise SD = 0, so there's
                                        // nothing to standardize
                    {
                        double sum1 = 0., sum2 = 0., sumsq1 = 0., sumsq2 = 0.;
                        for (const auto &p : thisbin)
                            {
                                if (isfinite(p.second.first))
                                    {
                                        sum1 += p.second.first;
                                        sumsq1 += pow(p.second.first, 2.);
                                    }
                                if (isfinite(p.second.second))
                                    {
                                        sum2 += p.second.second;
                                        sumsq2 += pow(p.second.second, 2.);
                                    }
                            }
                        double mean1
                            = sum1 / static_cast<double>(thisbin.size());
                        double mean2
                            = sum2 / static_cast<double>(thisbin.size());
                        double C = static_cast<double>(thisbin.size())
                                   / static_cast<double>(thisbin.size() - 1);
                        double var1
                            = C * (sumsq1 / static_cast<double>(thisbin.size())
                                   - pow(mean1, 2.));
                        double var2
                            = C * (sumsq2 / static_cast<double>(thisbin.size())
                                   - pow(mean2, 2.));
                        double sd1 = sqrt(var1), sd2 = sqrt(var2);
                        for (const auto &data : thisbin)
                            {
                                rv = maxabs(data.second.first, mean1, sd1, rv);
                                rv2 = maxabs(data.second.second, mean2, sd2,
                                             rv2);
                            }
                    }
            }
        if (ttlSNPs != __filtered.numsites())
            {
                throw runtime_error("Incorrect number of SNPs processed");
            }
        return make_pair(rv, rv2);
    }
}
