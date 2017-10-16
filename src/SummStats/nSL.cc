#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <algorithm>
#include <numeric>
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

    vector<tuple<double, double, uint32_t>>
    nSL_t(const SimData &d, const std::unordered_map<double, double> &gmap)
    {
        vector<size_t> core_snps(d.numsites());
        for (size_t core = 0; core < core_snps.size(); ++core)
            core_snps[core] = core;
        return nSL_t(d, core_snps, gmap);
    }

    vector<tuple<double, double, uint32_t>>
    nSL_t(const SimData &d, const std::vector<size_t> &core_snps,
          const std::unordered_map<double, double> &gmap)
    {
        vector<tuple<double, double, uint32_t>> rv(core_snps.size());
        using offset_type = SimData::const_site_iterator::difference_type;
        tbb::parallel_for(
            tbb::blocked_range<std::size_t>(0, rv.size()),
            [&rv, &d, &core_snps,
             &gmap](const tbb::blocked_range<std::size_t> &r) {
                for (std::size_t i = r.begin(); i < r.end(); ++i)
                    {
                        auto temp = __nSLdetails(core_snps[i], d, gmap);
                        offset_type offset = static_cast<offset_type>(i);
                        uint32_t dcount = static_cast<uint32_t>(
                            count((d.sbegin() + offset)->second.begin(),
                                  (d.sbegin() + offset)->second.end(), '1'));
                        rv[i]
                            = make_tuple(log(temp[0]) - log(temp[2]),
                                         log(temp[1]) - log(temp[3]), dcount);
                    }
            });
        return rv;
    }

    // double
    // update_return_value(vector<double> &binned_scores, const double
    // current_rv)
    //// normalize, etc.
    //{
    //    if (binned_scores.size() > 1)
    //        {
    //            double sum = accumulate(binned_scores.begin(),
    //                                    binned_scores.end(), 0.);
    //            double sumsq = accumulate(
    //                binned_scores.begin(), binned_scores.end(), 0.,
    //                [](double ss, double val) { return ss + pow(val, 2.0);
    //                });
    //            double mean = sum /
    //            static_cast<double>(binned_scores.size());
    //            double sd = sqrt(sumsq);
    //            transform(binned_scores.begin(), binned_scores.end(),
    //                      binned_scores.begin(), [mean, sd](double score) {
    //                          return (score - mean) / sd;
    //                      });
    //            auto me = max_element(
    //                binned_scores.begin(), binned_scores.end(),
    //                [](double a, double b) { return fabs(a) < fabs(b); });
    //            if (!isfinite(current_rv) || fabs(*me) > fabs(current_rv))
    //                {
    //                    return fabs(*me);
    //                }
    //        }
    //    return current_rv;
    //}
    /*
      Return max. abs value of standardized nSL and iHS, with the latter as
      defined by Ferrer-Admetella et al.
    */
    // pair<double, double>
    // snSL(const SimData &d, const double minfreq, const double binsize,
    //     bool filter_minor, const std::unordered_map<double, double> &gmap)
    //{
    //    if (d.empty())
    //        return make_pair(std::numeric_limits<double>::quiet_NaN(),
    //                         std::numeric_limits<double>::quiet_NaN());

    //    vector<unsigned> dcounts;
    //    dcounts.reserve(d.numsites());
    //    vector<size_t> core_snps;

    //    for (auto p = d.sbegin(); p != d.send(); ++p)
    //        {
    //            unsigned dcount = static_cast<unsigned>(
    //                count(p->second.begin(), p->second.end(), '1'));
    //            if (dcount && dcount < d.size())
    //                {
    //                    double f = static_cast<double>(dcount)
    //                               / static_cast<double>(d.size());
    //                    if (filter_minor)
    //                        {
    //                            f = min(f, 1. - f);
    //                            dcount = min(dcount,
    //                                         static_cast<unsigned>(d.size())
    //                                             - dcount);
    //                        }
    //                    if (f >= minfreq)
    //                        {
    //                            core_snps.push_back(
    //                                static_cast<size_t>(p - d.sbegin()));
    //                            dcounts.push_back(dcount);
    //                        }
    //                }
    //        }
    //    if (core_snps.empty())
    //        return make_pair(std::numeric_limits<double>::quiet_NaN(),
    //                         std::numeric_limits<double>::quiet_NaN());
    //    // Get the stats
    //    auto nSLstats = nSL_t(d, core_snps, gmap);
    //    const std::size_t nbins = static_cast<size_t>(
    //        std::ceil(((filter_minor) ? 0.5 : 1.0) / binsize));
    //    vector<vector<double>> binned_scores_nSL(nbins),
    //        binned_scores_iHS(nbins);
    //    double binsize_counts = (binsize * double(d.size()));
    //    for (size_t i = 0; i < core_snps.size(); ++i)
    //        {
    //            size_t bin
    //                = size_t(static_cast<double>(dcounts[i]) /
    //                binsize_counts);
    //            if (isfinite(get<0>(nSLstats[i])))
    //                {
    //                    binned_scores_nSL[bin].push_back(get<0>(nSLstats[i]));
    //                }
    //            if (isfinite(get<1>(nSLstats[i])))
    //                {
    //                    binned_scores_iHS[bin].push_back(get<1>(nSLstats[i]));
    //                }
    //        }
    //    double rv_nSL = numeric_limits<double>::quiet_NaN(),
    //           rv_iHS = numeric_limits<double>::quiet_NaN();
    //    for (size_t i = 0; i < binned_scores_nSL.size(); ++i)
    //        {
    //            rv_nSL = update_return_value(binned_scores_nSL[i], rv_nSL);
    //            rv_iHS = update_return_value(binned_scores_iHS[i], rv_iHS);
    //        }
    //    return make_pair(rv_nSL, rv_iHS);
    //}
}
