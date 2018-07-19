#include <cmath>
#include <limits>
#include <Sequence/summstats/thetapi.hpp>
#include <Sequence/summstats/thetaw.hpp>
#include <Sequence/summstats/nvariablesites.hpp>
#include <Sequence/summstats/auxillary.hpp>
namespace Sequence
{
    double
    tajd(const VariantMatrix& m)
    {
        auto S = total_number_of_mutations(m);
        if (S == 0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
        auto w = thetaw(m);
        auto pi = thetapi(m);
        auto a1 = summstats_aux::a_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto a2 = summstats_aux::b_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto dn = static_cast<double>(m.nsam);
        double b1 = (dn + 1.0) / (3.0 * (dn - 1.0));
        double b2
            = (2.0 * (std::pow(dn, 2.0) + dn + 3.0)) / (9.0 * dn * (dn - 1.0));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (dn + 2.0) / (a1 * dn) + a2 / std::pow(a1, 2.0);
        double e1 = c1 / a1;
        double e2 = c2 / (std::pow(a1, 2.0) + a2);
        double denominator = std::pow((e1 * S + e2 * S * (S - 1.0)), 0.5);
        return (pi - w) / denominator;
    }
} // namespace Sequence
