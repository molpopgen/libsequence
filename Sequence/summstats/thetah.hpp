#ifndef SEQUENCE_SUMMSTATS_THETAH_HPP__
#define SEQUENCE_SUMMSTATS_THETAH_HPP__

#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    double thetah(const VariantMatrix& m, const std::int8_t refstate);

    double thetah(const VariantMatrix& m,
                  const std::vector<std::int8_t>& refstates);
} // namespace Sequence

#endif
