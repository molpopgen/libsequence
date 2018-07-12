#include <unordered_map>
#include <stdexcept>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace Sequence
{
    double
    thetapi(const VariantMatrix& m)
    {
        double pi = 0.0;
        for (std::size_t site = 0; site < m.nsites; ++site)
            {
                auto site_view = get_RowView(m, site);
                std::unordered_map<std::int8_t, unsigned> counts;
                unsigned nmissing = 0;
                for (auto i : site_view)
                    {
                        if (i == VariantMatrix::mask)
                            {
                                throw std::runtime_error(
                                    "invalid state value encountered");
                            }
                        else if (i < 0)
                            {
                                ++nmissing;
                            }
                        else
                            {
                                counts[i]++;
                            }
                    }
                double sample_size
                    = static_cast<double>(site_view.row_size - nmissing);
                double sshom = 0.0; //sum of site homozygosity
                for (auto& x : counts)
                    {
                        sshom += static_cast<double>(x.second)
                                 * static_cast<double>(x.second - 1.0);
                    }
                pi += 1.0 - sshom / (sample_size * (sample_size - 1.0));
            }
        return pi;
    }

} // namespace Sequence
