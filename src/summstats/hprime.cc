#include <cmath>
#include <functional>
#include <Sequence/summstats/algorithm.hpp>
#include "hprime_faywuh_aggregator.hpp"
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace
{
    double
    hprime_common(const Sequence::VariantMatrix &m, const unsigned S,
                  const double tp, const double tl)
    {
        using namespace Sequence;
        if (tp == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto a = summstats_aux::a_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto b = summstats_aux::b_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto b1
            = summstats_aux::b_sub_n_plus1(static_cast<std::uint32_t>(m.nsam));

        double tw
            = static_cast<double>(S) / a; //TODO: replace with call to thetaw
        double tsq = S * (S - 1) / (a * a + b);
        double n = static_cast<double>(m.nsam);

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
    hprime(const VariantMatrix &m, const std::int8_t refstate)
    {
        detail::hprime_faywuh_aggregator hp(1.0);
        sstats_algo::aggregate_sites(m, std::ref(hp), refstate);
        return hprime_common(m, hp.S, hp.pi, hp.theta);
    }

    double
    hprime(const VariantMatrix &m, const std::vector<std::int8_t> &refstates)
    {
        detail::hprime_faywuh_aggregator hp(1.0);
        sstats_algo::aggregate_sites(m, std::ref(hp), refstates);
        return hprime_common(m, hp.S, hp.pi, hp.theta);
    }
} // namespace Sequence
