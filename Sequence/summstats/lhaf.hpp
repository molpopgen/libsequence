#ifndef SEQUENCE_SUMMSTATS_LHAF_HPP
#define SEQUENCE_SUMMSTATS_LHAF_HPP

#include <cstdint>
#include <vector>
#include <Sequence/VariantMatrix.hpp>

namespace Sequence
{
    /*! \brief l-Haf statistic of \cite Ronen2015-te
    * \param m A VariantMatrix
    * \param refstate The ancstral state
    * \param l The power parameter
    * \return vector of the statistic
    * \ingroup popgenanalysis
    */
    std::vector<double> lhaf(const VariantMatrix &m,
                             const std::int8_t refstate, const double l);
} // namespace Sequence
#endif
