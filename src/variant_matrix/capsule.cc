#include <Sequence/VectorCapsules.hpp>
#include <limits>
#include <algorithm>
#include <iterator>
#include <cmath>

namespace Sequence
{
    std::size_t&
    VectorGenotypeCapsule::nsites()
    {
        return nsites_;
    }

    std::size_t&
    VectorGenotypeCapsule::nsam()
    {
        return nsam_;
    }

    std::size_t
    VectorGenotypeCapsule::nsites() const
    {
        return nsites_;
    }

    std::size_t
    VectorGenotypeCapsule::nsam() const
    {
        return nsam_;
    }

    std::int8_t& VectorGenotypeCapsule::operator()(std::size_t site, std::size_t sample)
    {
        return buffer[site*nsam_ + sample];
    }

    const std::int8_t& VectorGenotypeCapsule::operator()(std::size_t site, std::size_t sample) const
    {
        return buffer[site*nsam_ + sample];
    }

    std::int8_t*
    VectorGenotypeCapsule::data()
    {
        return buffer.data();
    }

    const std::int8_t*
    VectorGenotypeCapsule::data() const
    {
        return buffer.data();
    }

    const std::int8_t*
    VectorGenotypeCapsule::cdata() const
    {
        return buffer.data();
    }

    std::unique_ptr<GenotypeCapsule>
    VectorGenotypeCapsule::clone() const
    {
        return std::unique_ptr<GenotypeCapsule>(
            new VectorGenotypeCapsule(this->buffer, this->nsites_));
    }

    std::int8_t*
    VectorGenotypeCapsule::begin()
    {
        return buffer.data();
    }

    const std::int8_t*
    VectorGenotypeCapsule::begin() const
    {
        return buffer.data();
    }

    std::int8_t*
    VectorGenotypeCapsule::end()
    {
        return buffer.data() + buffer.size();
    }

    const std::int8_t*
    VectorGenotypeCapsule::end() const
    {
        return buffer.data() + buffer.size();
    }

    const std::int8_t*
    VectorGenotypeCapsule::cbegin() const
    {
        return begin();
    }

    const std::int8_t*
    VectorGenotypeCapsule::cend() const
    {
        return end();
    }

    bool
    VectorGenotypeCapsule::empty() const
    {
        return buffer.empty();
    }

    std::size_t
    VectorGenotypeCapsule::size() const
    {
        return buffer.size();
    }

    bool
    VectorGenotypeCapsule::resizable() const
    {
        return true;
    }

    void
    VectorGenotypeCapsule::resize(bool remove_sites)
    {
        buffer.erase(std::remove_if(
                         std::begin(buffer), std::end(buffer),
                         [this](std::int8_t x) {
                             return x
                                    == std::numeric_limits<std::int8_t>::min();
                         }),
                     std::end(buffer));
        if (remove_sites)
            {
                nsites_ = buffer.size() / nsam_;
            }
        else
            {
                nsam_ = buffer.size() / nsites_;
            }
    }

    // Position capsule based on std::vector here
    double& VectorPositionCapsule::operator[](std::size_t i)
    {
        return buffer[i];
    }

    const double& VectorPositionCapsule::operator[](std::size_t i) const
    {
        return buffer[i];
    }

    double*
    VectorPositionCapsule::data()
    {
        return buffer.data();
    }

    const double*
    VectorPositionCapsule::data() const
    {
        return buffer.data();
    }

    const double*
    VectorPositionCapsule::cdata() const
    {
        return buffer.data();
    }

    std::unique_ptr<PositionCapsule>
    VectorPositionCapsule::clone() const
    {
        return std::unique_ptr<PositionCapsule>(
            new VectorPositionCapsule(*this));
    }

    double*
    VectorPositionCapsule::begin()
    {
        return buffer.data();
    }

    const double*
    VectorPositionCapsule::begin() const
    {
        return buffer.data();
    }

    double*
    VectorPositionCapsule::end()
    {
        return buffer.data() + buffer.size();
    }

    const double*
    VectorPositionCapsule::end() const
    {
        return buffer.data() + buffer.size();
    }

    const double*
    VectorPositionCapsule::cbegin() const
    {
        return begin();
    }

    const double*
    VectorPositionCapsule::cend() const
    {
        return end();
    }

    bool
    VectorPositionCapsule::empty() const
    {
        return buffer.empty();
    }

    std::size_t
    VectorPositionCapsule::size() const
    {
        return buffer.size();
    }

    std::size_t
    VectorPositionCapsule::nsites() const
    {
        return buffer.size();
    }

    bool
    VectorPositionCapsule::resizable() const
    {
        return true;
    }

    void
    VectorPositionCapsule::resize(bool /*unused*/)
    {
        buffer.erase(
            std::remove_if(std::begin(buffer), std::end(buffer),
                           [this](double d) { return std::isnan(d); }),
            std::end(buffer));
    }
} // namespace Sequence
