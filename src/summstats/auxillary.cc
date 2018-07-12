#include <cstdint>
#include <cmath>

namespace Sequence
{
    namespace summstats_aux
    {
        double
        a_sub_n(const std::uint32_t nsam)
        {
            double rv = 0.0;
            for (std::uint32_t i = 1; i < nsam; ++i)
                {
                    rv += 1.0 / static_cast<double>(i);
                }
            return rv;
        }
		
        double
        b_sub_n(const std::uint32_t nsam)
        {
            double rv = 0.0;
            for (std::uint32_t i = 1; i < nsam; ++i)
                {
                    rv += 1.0 / std::pow(static_cast<double>(i), 2.0);
                }
            return rv;
        }
    } // namespace summstats_aux
} // namespace Sequence
