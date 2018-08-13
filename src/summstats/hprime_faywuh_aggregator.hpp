#ifndef SEQUENCE_DETAIL_HPRIME_FAYWUH_AGGREGATOR_HPP
#define SEQUENCE_DETAIL_HPRIME_FAYWUH_AGGREGATOR_HPP

#include <Sequence/AlleleCountMatrix.hpp>
#include <stdexcept>
#include <cmath>

namespace Sequence
{
    namespace detail
    {
        struct stat_is_faywuh
        {
        };

        struct stat_is_hprime
        {
        };

        class hprime_faywuh_row_processor
        {
          private:
            inline double
            update_stat(const double d, const stat_is_faywuh)
            {
                return std::pow(d, 2.0);
            }

            inline double
            update_stat(const double d, const stat_is_hprime)
            {
                return d;
            }

            inline double
            denominator(const std::size_t nsam, const stat_is_faywuh)
            {
                return 2./static_cast<double>(nsam * (nsam - 1));
            }

            inline double
            denominator(const std::size_t nsam, const stat_is_hprime)
            {
                return 1./static_cast<double>(nsam - 1);
            }

            template <typename T>
            inline void
            stat_details(const Sequence::AlleleCountMatrix &ac,
                         const std::size_t i, const std::size_t refindex,
                         const T t)
            {
                unsigned nstates = 0;
                bool refseen = false;
                double temp = 0.0;
                double homozygosity = 0.0;
                for (std::size_t j = i; j < i + ac.ncol; ++j)
                    {
                        auto ci = ac.counts[j];
                        if (ci > 0)
                            {
                                homozygosity
                                    += static_cast<double>(ci * (ci - 1));
                                ++nstates;
                                if (j - i != refindex)
                                    {
                                        temp += update_stat(ci, t);
                                    }
                                else
                                    {
                                        refseen = true;
                                    }
                            }
                    }
                if (nstates > 2)
                    {
                        throw std::runtime_error(
                            "site has more than one derived state");
                    }

                if (nstates > 1)
                    {
                        ++S;
                    }
                if (refseen)
                    {
                        double nnm1
                            = static_cast<double>(ac.nsam * (ac.nsam - 1));
                        pi += 1.0 - homozygosity / nnm1;
                        double x = denominator(ac.nsam, t);
                        theta += temp * x;
                    }
            }

          public:
            unsigned S;
            double pi, theta;
            hprime_faywuh_row_processor()
                : S{ 0 }, pi{ 0.0 }, theta{ 0.0 }
            {
            }

            inline void
            operator()(const Sequence::AlleleCountMatrix &ac,
                       const std::size_t i, const std::size_t refindex,
                       const stat_is_faywuh dispatch)
            {
                stat_details(ac, i, refindex, dispatch);
            }

            inline void
            operator()(const Sequence::AlleleCountMatrix &ac,
                       const std::size_t i, const std::size_t refindex,
                       const stat_is_hprime dispatch)
            {
                stat_details(ac, i, refindex, dispatch);
            }
        };
    } // namespace detail
} // namespace Sequence

#endif
