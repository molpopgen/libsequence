#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/thetapi.hpp>
#include <Sequence/summstats/thetal.hpp>
#include <Sequence/summstats/allele_counts.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    hprime(const VariantMatrix &m, const std::int8_t refstate)
    {
        auto tp = thetapi(m, refstate);
        auto tl = thetal(m, refstate);

        auto a = summstats_aux::a_sub_n(m.nsam);
        auto b = summstats_aux::b_sub_n(m.nsam);
        auto b1 = summstats_aux::b_sub_n_plus1(m.nsam);

        //Count number of bi-allelic sites
        //TODO: generalize this
    }
} // namespace Sequence
