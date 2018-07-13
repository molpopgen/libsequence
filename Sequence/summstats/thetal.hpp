#ifndef SEQUENCE_SUMMSTATS_THETAH_HPP__
#define SEQUENCE_SUMMSTATS_THETAH_HPP__

#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    double thetal(const VariantMatrix& m, const std::int8_t refstate);

    double thetal(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);
} // namespace Sequence

#endif
