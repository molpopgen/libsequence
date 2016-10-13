#include <Sequence/SimData.hpp>
#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace std;

namespace
{
    double
    update_s2(std::string::const_iterator start, std::string::const_iterator left,
              std::string::const_iterator right, const Sequence::SimData &d,
              const std::unordered_map<double, double> &gmap)
    {
        auto p1 = d.position(
            std::vector<double>::size_type(distance(start, right)) - 1);
        auto p2
            = d.position(std::vector<double>::size_type(distance(start, left)));
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
    pair<double, double>
    __nlSsum(const unsigned &core, const Sequence::SimData &d,
             const vector<size_t> &coretype,
             const std::unordered_map<double, double> &gmap)
    {
        double s = 0., s2 = 0.;
        unsigned nc = 0u;
        auto csize = coretype.size();
        for (size_t i = 0; i < csize; ++i)
            {
				auto start = d[coretype[i]].cbegin();
                auto bi =  start + core;
                auto eri = d[coretype[i]].crend();
                auto ei = d[coretype[i]].cend();
                for (size_t j = i + 1; j < csize; ++j)
                    {
                        auto bj = d[coretype[j]].cbegin() + core;
                        auto right = mismatch(bi, ei, bj);
                        string::const_reverse_iterator ri1(bi), ri2(bj);
                        auto left = mismatch(ri1, eri, ri2);
                        if (left.first != eri && right.first != ei)
                            {
                                s += static_cast<double>(
                                    distance(left.first.base(), right.first)
                                    + 1);
                                s2 += update_s2(start, left.first.base(), right.first,
                                                d, gmap);
                                ++nc;
                            }
                    }
            }
        return make_pair(s / static_cast<double>(nc),
                         s2 / static_cast<double>(nc));
    }
}

namespace Sequence
{
    /*
      The nSL statistic of doi:10.1093/molbev/msu077
    */
    pair<double, double>
    nSL(const unsigned &core, const SimData &d,
        const std::unordered_map<double, double> &gmap
        = std::unordered_map<double, double>())
    {
        std::vector<size_t> der, anc;
        der.reserve(d.size());
        anc.reserve(d.size());
        for (unsigned i = 0; i < d.size(); ++i)
            {
                if (d[i][core] == '1')
                    der.push_back(i);
                else
                    anc.push_back(i);
            }
        pair<double, double> A = __nlSsum(core, d, anc, gmap),
                             D = __nlSsum(core, d, der, gmap);
        return make_pair(log(A.first) - log(D.first),
                         log(A.second) - log(D.second));
    }

    /*
      Return max. abs value of standardized nSL and iHS, with the latter as
      defined by Ferrer-Admetella et al.
    */
    pair<double, double>
    snSL(const SimData &d, const double minfreq, const double binsize,
         const std::unordered_map<double, double> &gmap
         = std::unordered_map<double, double>())
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
        // Associate the stats with their DAFs
        vector<pair<double, pair<double, double>>> binning;
        for (unsigned i = 0; i < __filtered.numsites(); ++i)
            {
                pair<double, double> rvi = nSL(i, __filtered, gmap);
                binning.push_back(
                    make_pair(static_cast<double>(dcounts[i])
                                  / static_cast<double>(__filtered.size()),
                              rvi));
            }
        double rv = std::numeric_limits<double>::quiet_NaN(),
               rv2 = std::numeric_limits<double>::quiet_NaN();
        // Now, bin, standardise, and move on...
        for (double l = minfreq; l < 1.; l += binsize)
            {
                vector<pair<double, pair<double, double>>> thisbin;
                copy_if(binning.begin(), binning.end(), back_inserter(thisbin),
                        [&](const pair<double, pair<double, double>> &data) {
                            return isfinite(data.second.first)
                                   && data.first >= l
                                   && data.first < l + binsize;
                        });
                if (thisbin.size() > 1) // otherwise SD = 0, so there's
                                        // nothing to standardize
                    {
                        double mean1 = 0., mean2 = 0.;
                        for (const auto &p : thisbin)
                            {
                                if (isfinite(p.second.first))
                                    mean1 += p.second.first;
                                if (isfinite(p.second.second))
                                    mean2 += p.second.second;
                            }
                        mean1 /= static_cast<double>(thisbin.size());
                        mean2 /= static_cast<double>(thisbin.size());
                        double var1 = 0., var2 = 0.;
                        for (const auto &p : thisbin)
                            {
                                if (isfinite(p.second.first))
                                    var1 += pow(p.second.first - mean1, 2.0);
                                if (isfinite(p.second.first))
                                    var2 += pow(p.second.second - mean2, 2.0);
                            }
                        var1 /= static_cast<double>(thisbin.size() - 1);
                        var2 /= static_cast<double>(thisbin.size() - 1);
                        double sd1 = sqrt(var1), sd2 = sqrt(var2);
                        for_each(
                            thisbin.begin(), thisbin.end(),
                            [&](const pair<double, pair<double, double>>
                                    &data) {
                                double z1
                                    = (isfinite(sd1))
                                          ? (data.second.first - mean1) / sd1
                                          : numeric_limits<double>::
                                                quiet_NaN(),
                                    z2
                                    = (isfinite(sd2))
                                          ? (data.second.second - mean2) / sd2
                                          : numeric_limits<double>::
                                                quiet_NaN();
                                // Get max abs val of each stat
                                if (isfinite(z1)
                                    && (!isfinite(rv) || fabs(z1) > fabs(rv)))
                                    rv = z1;
                                if (isfinite(z2) && (!isfinite(rv2)
                                                     || fabs(z2) > fabs(rv2)))
                                    rv2 = z2;
                            });
                    }
            }
        return make_pair(rv, rv2);
    }
}
