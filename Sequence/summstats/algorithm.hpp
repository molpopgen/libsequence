/*!
\file summstats/algorithm.hpp
*/
#ifndef SEQUENCE_SUMMSTATS_ALGORITHM_HPP
#define SEQUENCE_SUMMSTATS_ALGORITHM_HPP

#include <cstdint>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/StateCounts.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    namespace sstats_algo
    {
        template <typename F>
        inline void
        aggregate_sites(const VariantMatrix& m, const F& f,
                        const std::int8_t refstate)
        /*! Helper algorithm for implementing summary statistics.
         *
         *  Several common summary statistics are combinations of
         *  others.  Examples include Tajima's D, Fay and Wu's H,
         *  etc..  If we take D as an example, it is tempting to
         *  use existing functions, such as Sequence::thetapi and
         *  Sequence::thetaw, as intermediate steps.
         *  However, doing so goes over the data multiple times.  
         *
         *  Fortunately, these statistics are often easy enough to 
         *  implement that we could calculate pi and Watterson's theta 
         *  in one loop.  This function helps you do that.
         *  
         *  \param m A VariantMatrix
         *  \param f A function taking a const StateCounts & and returning nothing.
         *  \param refstate The reference state.
         *
         *  This function loops over \a m.nsites and passes the state counts on
         *  to the aggregator function \a f.
         *
         *  See the implementation of Sequence::tajd for an example.
         *
         *  \ingroup popgenanalysis
         */
        {
            StateCounts c(refstate);
            for (std::size_t i = 0; i < m.nsites; ++i)
                {
                    auto r = get_ConstRowView(m, i);
                    c(r);
                    f(c);
                }
        }
    } // namespace sstats_algo
} // namespace Sequence

#endif
