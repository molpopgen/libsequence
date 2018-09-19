#ifndef SEQUENCE_SUMMSTATS_ALGORITHM
#define SEQUENCE_SUMMSTATS_ALGORITHM

#include <cstdint>
#include <utility>
#include <algorithm>

namespace Sequence
{
    namespace summstats_algo
    {
        template <typename iterator>
        inline std::pair<iterator, iterator>
        mismatch_skip_missing(iterator beg, iterator end, iterator beg2)
        {
            auto m = std::mismatch(beg, end, beg2);
            while (m.first < end && (*m.first < 0 || *m.second < 0))
                {
                    m = std::mismatch(m.first + 1, end, m.second + 1);
                }
            return m;
        }

        template <typename iterator>
        inline std::int32_t
        ndiff_skip_missing(iterator beg, iterator end, iterator beg2)
        {
            std::int32_t ndiffs = 0;
            auto m = mismatch_skip_missing(beg, end, beg2);

            while (m.first < end)
                {
                    ++ndiffs;
                    m = mismatch_skip_missing(m.first + 1, end, m.second + 1);
                }
            return ndiffs;
        }
    } // namespace summstats_algo
} // namespace Sequence

#endif
