#include <cstdint>
#include <Sequence/summstats/ld.hpp>
#include <Sequence/summstats/allele_counts.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    std::int32_t
    rmin(const VariantMatrix& m)
    {
        if (m.nsites < 2)
            {
                return -1;
            }
        auto ac = allele_counts(m);
        std::vector<std::size_t> biallelic_site_indexes;
        for (std::size_t i = 0; i < ac.size(); ++i)
            {
                if (ac[i].nstates == 2)
                    {
                        biallelic_site_indexes.push_back(i);
                    }
            }
        if (biallelic_site_indexes.size() < 2)
            {
                return -1;
            }
        bool flag = false;
        std::size_t x = 0;
        std::int32_t rv = 0;
        for (std::size_t a = x + 1; a < biallelic_site_indexes.size(); ++a)
            {
                for (std::size_t b = (!flag) ? x : a - 1; b < a; ++b)
                    {
                        flag = false;
                        // We do not allow missing data to result in
                        // additional haplotypes
                        auto tl = two_locus_haplotype_counts(m, a, b, true);
                        if (tl.size() == 4)
                            {
                                ++rv;
                                flag = true;
                                break;
                            }
                    }
                if (flag == true)
                    {
                        x = a;
                    }
            }
        return rv;
    }
} // namespace Sequence
