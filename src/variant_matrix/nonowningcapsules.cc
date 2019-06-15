#include <Sequence/NonOwningCapsules.hpp>
#include <limits>
#include <algorithm>
#include <iterator>
#include <exception>

namespace
{
    void
    raise()
    {
        throw std::runtime_error("data are read-only");
    }
} // namespace

namespace Sequence
{
    std::size_t&
    NonOwningGenotypeCapsule::nsites()
    {
        return nsites_;
    }

    std::size_t&
    NonOwningGenotypeCapsule::nsam()
    {
        return nsam_;
    }

    std::size_t
    NonOwningGenotypeCapsule::nsites() const
    {
        return nsites_;
    }

    std::size_t
    NonOwningGenotypeCapsule::nsam() const
    {
        return nsam_;
    }

    std::size_t
    NonOwningGenotypeCapsule::row_offset() const
    {
        return k1;
    }

    std::size_t
    NonOwningGenotypeCapsule::col_offset() const
    {
        return k2;
    }

    std::int8_t&
    NonOwningGenotypeCapsule::operator()(std::size_t site, std::size_t sample)
    {
        raise();
        return *const_cast<std::int8_t*>(
            &buffer[(k1 * tda + k2) + site * tda + sample]);
    }

    const std::int8_t&
    NonOwningGenotypeCapsule::operator()(std::size_t site,
                                         std::size_t sample) const
    {
        //m'(i,j) = m->data[(k1*m->tda + k2) + i*m->tda + j]
        return buffer[(k1 * tda + k2) + site * tda + sample];
    }

    std::int8_t*
    NonOwningGenotypeCapsule::data()
    {
        raise();
        return const_cast<std::int8_t*>(buffer.get());
    }

    const std::int8_t*
    NonOwningGenotypeCapsule::data() const
    {
        return buffer.get();
    }

    const std::int8_t*
    NonOwningGenotypeCapsule::cdata() const
    {
        return buffer.get();
    }

    std::unique_ptr<GenotypeCapsule>
    NonOwningGenotypeCapsule::clone() const
    {
        return std::unique_ptr<GenotypeCapsule>(new NonOwningGenotypeCapsule(
            this->buffer.get(), this->nsites_, this->nsam_, this->k1, this->k2,
            this->tda));
    }

    std::int8_t*
    NonOwningGenotypeCapsule::begin()
    {
        raise();
        return const_cast<std::int8_t*>(buffer.get());
    }

    const std::int8_t*
    NonOwningGenotypeCapsule::begin() const
    {
        return buffer.get();
    }

    std::int8_t*
    NonOwningGenotypeCapsule::end()
    {
        raise();
        return const_cast<std::int8_t*>(buffer.get()) + nsam_ * nsites_;
    }

    const std::int8_t*
    NonOwningGenotypeCapsule::end() const
    {
        return buffer.get() + nsam_ * nsites_;
    }

    const std::int8_t*
    NonOwningGenotypeCapsule::cbegin() const
    {
        return cdata();
    }

    const std::int8_t*
    NonOwningGenotypeCapsule::cend() const
    {
        return cdata() + nsam_ + nsites_;
    }

    bool
    NonOwningGenotypeCapsule::empty() const
    {
        return nsam_ == 0 || nsites_ == 0;
    }

    std::size_t
    NonOwningGenotypeCapsule::size() const
    {
        return nsam_ * nsites_;
    }

    bool
    NonOwningGenotypeCapsule::resizable() const
    {
        return false;
    }

    // Position capsule based on std::vector here
    double& NonOwningPositionCapsule::operator[](std::size_t i)
    {
        raise();
        return *const_cast<double*>(&buffer[i]);
    }

    const double& NonOwningPositionCapsule::operator[](std::size_t i) const
    {
        return buffer[i];
    }

    double*
    NonOwningPositionCapsule::data()
    {
        raise();
        return const_cast<double*>(buffer.get());
    }

    const double*
    NonOwningPositionCapsule::data() const
    {
        return buffer.get();
    }

    const double*
    NonOwningPositionCapsule::cdata() const
    {
        return buffer.get();
    }

    std::unique_ptr<PositionCapsule>
    NonOwningPositionCapsule::clone() const
    {
        return std::unique_ptr<PositionCapsule>(
            new NonOwningPositionCapsule(buffer.get(), current_size));
    }

    double*
    NonOwningPositionCapsule::begin()
    {
        raise();
        return const_cast<double*>(buffer.get());
    }

    const double*
    NonOwningPositionCapsule::begin() const
    {
        return buffer.get();
    }

    double*
    NonOwningPositionCapsule::end()
    {
        raise();
        return const_cast<double*>(cdata()) + current_size;
    }

    const double*
    NonOwningPositionCapsule::end() const
    {
        return buffer.get() + current_size;
    }

    const double*
    NonOwningPositionCapsule::cbegin() const
    {
        return cdata();
    }

    const double*
    NonOwningPositionCapsule::cend() const
    {
        return cdata() + current_size;
    }

    bool
    NonOwningPositionCapsule::empty() const
    {
        return current_size == 0;
    }

    std::size_t
    NonOwningPositionCapsule::size() const
    {
        return current_size;
    }

    std::size_t
    NonOwningPositionCapsule::nsites() const
    {
        return current_size;
    }

    bool
    NonOwningPositionCapsule::resizable() const
    {
        return false;
    }
} // namespace Sequence

