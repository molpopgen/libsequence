#include <cstdint>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrix.hpp>

namespace
{
    inline double
    update(const Sequence::StateCounts& c, const std::int8_t refstate,
           const double power)
    {
        double H = 0.0;
        double nnm1 = static_cast<double>(c.n * (c.n - 1));
        std::size_t has_missing = (c.counts.find(-1) != c.counts.end());
        std::size_t num_non_missing_states = c.counts.size() - has_missing;
        bool ref_seen = false;
        for (auto& x : c.counts)
            {
                if (!(x.first < 0))
                    {
                        if (x.first != refstate)
                            {
                                H += std::pow(x.second, power);
                            }
                        else
                            {
                                ref_seen = true;
                            }
                    }
            }
        if (!ref_seen && num_non_missing_states > 1)
            {
                throw std::runtime_error(
                    "the site has more than one derived state");
            }
        if (!ref_seen) //The reference state must still be segregating
                       // TODO: need unit test for this case
            {
                return 0.0;
            }
        H /= nnm1;
        return H;
    }

    inline double
    calc_stat(const Sequence::VariantMatrix& m, const std::int8_t refstate,
              const double power)
    {
        if (refstate < 0)
            {
                throw std::invalid_argument(
                    "all reference states are encoded as missing");
            }
        double stat = 0.0;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                Sequence::StateCounts counts(r, refstate);
                stat += update(counts, refstate, power);
            }
        return stat;
    }

    inline double
    calc_stat(const Sequence::VariantMatrix& m,
              const std::vector<std::int8_t>& refstates, const double power)
    {
        if (std::all_of(refstates.begin(), refstates.end(),
                        [](const std::int8_t x) { return x < 0; }))
            {
                throw std::invalid_argument(
                    "all reference states are encoded as missing");
            }
        double stat = 0.0;
        for (std::size_t i = 0; i < m.nsites; ++i)
            {
                auto r = get_ConstRowView(m, i);
                Sequence::StateCounts counts(r, refstates[i]);
                stat += update(counts, refstates[i], power);
            }
        return stat;
    }
} // namespace

namespace Sequence
{
    double
    thetah(const VariantMatrix& m, const std::int8_t refstate)
    {
        return calc_stat(m, refstate, 2.0);
    }

    double
    thetah(const VariantMatrix& m, const std::vector<std::int8_t>& refstates)
    {
        return calc_stat(m, refstates, 2.0);
    }

    double
    thetal(const VariantMatrix& m, const std::int8_t refstate)
    {
        return calc_stat(m, refstate, 1.0);
    }

    double
    thetal(const VariantMatrix& m, const std::vector<std::int8_t>& refstates)
    {
        return calc_stat(m, refstates, 1.0);
    }
} // namespace Sequence
