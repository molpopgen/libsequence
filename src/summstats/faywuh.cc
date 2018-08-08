#include <functional>
#include <Sequence/summstats/algorithm.hpp>
#include "hprime_faywuh_aggregator.hpp"

namespace Sequence
{
    double
    faywuh(const VariantMatrix& m, const std::int8_t refstate)
    {
        detail::hprime_faywuh_aggregator agg(2.0);
        sstats_algo::aggregate_sites(m, std::ref(agg), refstate);
        if (agg.pi == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        return agg.pi - agg.theta;
    }

    double
    faywuh(const VariantMatrix& m, const std::vector<std::int8_t>& refstates)
    {
        detail::hprime_faywuh_aggregator agg(2.0);
        sstats_algo::aggregate_sites(m, std::ref(agg), refstates);
        if (agg.pi == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        return agg.pi - agg.theta;
    }
    double
    faywuh(const AlleleCountMatrix& ac, const std::int8_t refstate)
    {
        auto refindex = static_cast<std::size_t>(refstate);
        if (refindex >= ac.row_size)
            {
                throw std::invalid_argument(
                    "reference state greater than max allelic state");
            }
        detail::hprime_faywuh_row_processor rp(2.0);
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.row_size)
            {
                rp(ac, i, refindex);
            }
        if (!rp.S)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        return rp.pi - rp.theta;
    }

    double
    faywuh(const AlleleCountMatrix& ac,
           const std::vector<std::int8_t>& refstates)
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

        detail::hprime_faywuh_row_processor rp(2.0);
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
        if (!rp.S)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        return rp.pi - rp.theta;
    }
} // namespace Sequence
