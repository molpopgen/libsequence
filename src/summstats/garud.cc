#include <Sequence/summstats/garud.hpp>

namespace Sequence
{
    GarudStats::GarudStats()
        : H1(1.), H12(std::numeric_limits<double>::quiet_NaN()),
          H2H1(std::numeric_limits<double>::quiet_NaN())
    {
    }

    GarudStats::GarudStats(const double __h1, const double __h12,
                           const double __h2h1)
        : H1(__h1), H12(__h12), H2H1(__h2h1)
    {
    }
} // namespace Sequence
