#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/thetapi.hpp>
#include <Sequence/summstats/thetal.hpp>
#include <Sequence/summstats/nvariablesites.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    hprime(const VariantMatrix &m, const std::int8_t refstate)
    {
        auto tp = thetapi(m);
        auto tl = thetal(m, refstate);
        if (tp == 0.0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

        auto a = summstats_aux::a_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto b = summstats_aux::b_sub_n(static_cast<std::uint32_t>(m.nsam));
        auto b1
            = summstats_aux::b_sub_n_plus1(static_cast<std::uint32_t>(m.nsam));

        //Count number of bi-allelic sites
        //TODO: generalize this
        auto S = nbiallelic_sites(m);
    }
} // namespace Sequence
