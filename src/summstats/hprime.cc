#include <cmath>
#include <functional>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/summstats/auxillary.hpp>
#include "hprime_faywuh_aggregator.hpp"

namespace
{
    double
    hprime_common(const std::uint32_t nsam, const unsigned S, const double tp,
                  const double tl)
    {
        using namespace Sequence;
        if (tp == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto a = summstats_aux::a_sub_n(static_cast<std::uint32_t>(nsam));
        auto b = summstats_aux::b_sub_n(static_cast<std::uint32_t>(nsam));
        auto b1
            = summstats_aux::b_sub_n_plus1(static_cast<std::uint32_t>(nsam));

        double tw
            = static_cast<double>(S) / a; //TODO: replace with call to thetaw
        double tsq = S * (S - 1) / (a * a + b);
        double n = static_cast<double>(nsam);

        double vThetal
            = (n * tw) / (2.0 * (n - 1.0))
              + (2.0 * std::pow(n / (n - 1.0), 2.0) * (b1 - 1.0) - 1.0) * tsq;
        double vPi = (3.0 * n * (n + 1.0) * tw + 2.0 * (n * n + n + 3.0) * tsq)
                     / (9 * n * (n - 1.0));
        double cov
            = ((n + 1.0) / (3.0 * (n - 1.0))) * tw
              + ((7.0 * n * n + 3.0 * n - 2.0 - 4.0 * n * (n + 1.0) * b1)
                 / (2.0 * std::pow((n - 1.0), 2.0)))
                    * tsq;
        return (tp - tl) / std::pow(vThetal + vPi - 2.0 * cov, 0.5);
    }
} // namespace

namespace Sequence
{
    double
    hprime(const AlleleCountMatrix &ac, const std::int8_t refstate)
    {
        auto refindex = static_cast<std::size_t>(refstate);
        if (refindex >= ac.row_size)
            {
                throw std::invalid_argument(
                    "reference state greater than max allelic state");
            }
        detail::hprime_faywuh_row_processor rp(1.0);
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.row_size)
            {
                rp(ac, i, refindex);
            }
        return hprime_common(ac.nsam, rp.S, rp.pi, rp.theta);
    }

    double
    hprime(const AlleleCountMatrix &ac,
           const std::vector<std::int8_t> &refstates)
    {
        if (ac.counts.empty())
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        if (refstates.size() != ac.counts.size() / ac.row_size)
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

        detail::hprime_faywuh_row_processor rp(1.0);
        std::size_t rstate = 0;
        for (std::size_t i = 0; i < ac.counts.size();
             i += ac.row_size, ++rstate)
            {
                if (refstates[rstate] >= 0)
                    {
                        auto refindex
                            = static_cast<std::size_t>(refstates[rstate]);
                        if (refindex >= ac.row_size)
                            {
                                throw std::invalid_argument(
                                    "reference state greater than max allelic "
                                    "state");
                            }
                        rp(ac, i, refindex);
                    }
            }
        return hprime_common(ac.nsam, rp.S, rp.pi, rp.theta);
    }
} // namespace Sequence
