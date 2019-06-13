#include <Sequence/VectorCapsules.hpp>

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

    std::int8_t& VectorGenotypeCapsule::operator[](std::size_t i)
    {
        return buffer[i];
    }

    const std::int8_t& VectorGenotypeCapsule::operator[](std::size_t i) const
    {
        return buffer[i];
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
} // namespace Sequence
