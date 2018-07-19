#ifndef SEQUENCE_SUMMSTATS_AUXILLARY_HPP__
#define SEQUENCE_SUMMSTATS_AUXILLARY_HPP__

#include <cstdint>

namespace Sequence
{
    namespace summstats_aux
    {
        double a_sub_n(const std::uint32_t);
        double b_sub_n(const std::uint32_t nsam);
        double b_sub_n_plus1(const std::uint32_t nsam);
    } // namespace summstats_aux
} // namespace Sequence

#endif
