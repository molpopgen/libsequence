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
} // namespace Sequence
