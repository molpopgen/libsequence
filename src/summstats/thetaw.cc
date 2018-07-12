#include <unordered_map>
#include <stdexcept>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/summstats/auxillary.hpp>

namespace Sequence
{
    double
    thetaw(const VariantMatrix& m)
    {
        double w = 0.0;
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
                if (counts.size() > 1)
                    {
                        auto nstates = counts.size();
                        std::uint32_t sample_size = static_cast<std::uint32_t>(
                            site_view.row_size - nmissing);
                        auto denom = summstats_aux::a_sub_n(sample_size);
						w += static_cast<double>(nstates-1)/denom;
                    }
            }
        return w;
    }

} // namespace Sequence
