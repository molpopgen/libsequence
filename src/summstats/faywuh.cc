#include <Sequence/summstats/thetah.hpp>
#include <Sequence/summstats/thetapi.hpp>

namespace Sequence
{
    double
    faywuh(const VariantMatrix& m, const std::int8_t refstate)
    {
        auto pi = thetapi(m);
        auto h = thetah(m, refstate);
        if (pi == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        return pi - h;
    }

    double
    faywuh(const VariantMatrix& m, const std::vector<std::int8_t>& refstates)
    {
        auto pi = thetapi(m);
        auto h = thetah(m, refstates);
        if (pi == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        return pi - h;
    }
} // namespace Sequence
