#ifndef SEQUENCE_DETAIL_HPRIME_FAYWUH_AGGREGATOR_HPP
#define SEQUENCE_DETAIL_HPRIME_FAYWUH_AGGREGATOR_HPP

#include <cmath>

namespace Sequence
{
    namespace detail
    {
        struct hprime_faywuh_row_processor
        {
            unsigned S;
            double pi, theta;
            const double power;
            hprime_faywuh_row_processor(const double p)
                : S{ 0 }, pi{ 0.0 }, theta{ 0.0 }, power{ p }
            {
            }

            inline void
            operator()(const Sequence::AlleleCountMatrix &ac,
                       const std::size_t i, const std::size_t refindex)
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
                                        temp += std::pow(
                                            static_cast<double>(ci), power);
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
                        double x = (power==1.0) ? 1.0/static_cast<double>(ac.nsam) : 
                            2.0/static_cast<double>(ac.nsam*(ac.nsam-1));
                        theta += temp * x;
                    }
            }
        };
    } // namespace detail
} // namespace Sequence

#endif
