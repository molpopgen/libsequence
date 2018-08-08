#include <stdexcept>
#include <Sequence/AlleleCountMatrix.hpp>

namespace Sequence
{
    double
    thetapi(const AlleleCountMatrix& ac)
    {
        double pi = 0.0;
        for (std::size_t i = 0; i < ac.counts.size(); i += ac.ncol)
            {
                std::int32_t nsam = 0;
                double homozygosity = 0.0;
                for (std::size_t j = i; j < i + ac.ncol; ++j)
                    {
                        nsam += ac.counts[j];
                        homozygosity += static_cast<double>(
                            ac.counts[j] * (ac.counts[j] - 1));
                    }
                pi += 1.0
                      - homozygosity / static_cast<double>(nsam * (nsam - 1));
            }
        return pi;
    }

} // namespace Sequence
