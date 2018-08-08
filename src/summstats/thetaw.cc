#include <algorithm>
#include <stdexcept>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    thetaw(const AlleleCountMatrix& ac)
    {
        double w = 0.0;
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.row_size)
            {
                std::int32_t nsam = 0, nstates = 0;
                for (std::size_t j = i; j < i + ac.row_size; ++j)
                    {
                        if (ac.counts[j] > 0)
                            {
                                nsam += ac.counts[j];
                                nstates++;
                            }
                    }
                if (nstates > 1)
                    {
                        auto denom = summstats_aux::a_sub_n(nsam);
                        w += static_cast<double>(nstates - 1) / denom;
                    }
            }
        return w;
    }

} // namespace Sequence
