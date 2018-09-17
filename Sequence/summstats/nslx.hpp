/// \file Sequence/summstats/nslx.hpp
/// \brief nSL and iHS
#ifndef SEQUENCE_SUMMSTATS_NSLX_HPP
#define SEQUENCE_SUMMSTATS_NSLX_HPP

#include <vector>
#include <cstdint>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    struct nSLiHS; //Forward declaration. See Sequence/summstats/nsl.hpp

    std::vector<nSLiHS> nslx(const VariantMatrix& m,
                             const std::int8_t refstate, const int x);
} // namespace Sequence

#endif
