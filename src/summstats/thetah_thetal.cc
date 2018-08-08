#include <cstdint>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <Sequence/AlleleCountMatrix.hpp>

namespace
{
    inline double
    process_row(const Sequence::AlleleCountMatrix& ac, const std::size_t i,
                const std::size_t refindex, const double power)
    {
        std::int32_t nsam = 0;
        std::int32_t nnonref = 0;
        double temp = 0.0;
        bool ref_seen = false;
        for (std::size_t j = i; j < i + ac.ncol; ++j)
            {
                if (ac.counts[j])
                    {
                        nsam += ac.counts[j];
                        if (j - i != refindex)
                            {
                                ++nnonref;
                                temp += std::pow(ac.counts[j], power);
                            }
                        else
                            {
                                ref_seen = true;
                            }
                    }
            }
        if (nnonref > 1)
            {
                throw std::runtime_error(
                    "site has more than one derived state");
            }
        if (ref_seen)
            {
                double nnm1 = static_cast<double>(nsam * (nsam - 1));
                return temp / nnm1;
            }
        return 0.0;
    }

    inline double
    calc_stat(const Sequence::AlleleCountMatrix& ac,
              const std::vector<std::int8_t>& refstates, const double power)
    {
        double rv = 0.0;
        if (ac.counts.empty())
            {
                return rv;
            }
        if (refstates.size() != ac.counts.size() / ac.ncol)
            {
                throw std::invalid_argument(
                    "incorrect number of reference states");
            }
        if (std::all_of(refstates.begin(), refstates.end(),
                        [](const std::int8_t i) { return i < 0; }))
            {
                throw std::invalid_argument(
                    "all reference states encoded as missing");
            }
        std::size_t rstate = 0;
        for (std::size_t i = 0; i < ac.counts.size();
             i += ac.ncol, ++rstate)
            {
                if (refstates[rstate] >= 0)
                    {
                        auto refindex
                            = static_cast<std::size_t>(refstates[rstate]);
                        if (refindex >= ac.ncol)
                            {
                                throw std::invalid_argument(
                                    "reference state greater than max allelic "
                                    "state");
                            }
                        rv += process_row(ac, i, refindex, power);
                    }
            }
        return rv;
    }

    inline double
    calc_stat(const Sequence::AlleleCountMatrix& ac,
              const std::int8_t refstate, const double power)
    {
        double rv = 0.0;
        if (ac.counts.empty())
            {
                return rv;
            }
        auto refindex = static_cast<std::size_t>(refstate);
        if (refindex >= ac.ncol)
            {
                throw std::invalid_argument(
                    "reference state greater than max allelic state");
            }
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.ncol)
            {
                rv += process_row(ac, i, refindex, power);
            }
        return rv;
    }
} // namespace

namespace Sequence
{
    double
    thetah(const AlleleCountMatrix& ac, const std::int8_t refstate)
    {
        return calc_stat(ac, refstate, 2.0);
    }

    double
    thetah(const AlleleCountMatrix& ac,
           const std::vector<std::int8_t>& refstates)
    {
        return calc_stat(ac, refstates, 2.0);
    }

    double
    thetal(const AlleleCountMatrix& ac, const std::int8_t refstate)
    {
        return calc_stat(ac, refstate, 1.0);
    }

    double
    thetal(const AlleleCountMatrix& ac,
           const std::vector<std::int8_t>& refstates)
    {
        return calc_stat(ac, refstates, 1.0);
    }

} // namespace Sequence
