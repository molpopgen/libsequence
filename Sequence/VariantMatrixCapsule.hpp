#ifndef SEQUENCE_VARIANT_MATRIX_CAPSULE
#define SEQUENCE_VARIANT_MATRIX_CAPSULE

#include <cstdint>
#include <memory>
#include <vector>

namespace Sequence
{
    template <typename T> struct Capsule
    {
        virtual ~Capsule() = default;
        // Following two may not be needed
        //virtual T& get(std::size_t site, std::size_t sample) = 0;
        //virtual const T& get(std::size_t site,
        //                               std::size_t sample) const = 0;
        virtual T* data() = 0;
        virtual const T* data() const = 0;
        virtual const T* cdata() const = 0;
        virtual T* begin() = 0;
        virtual const T* begin() const = 0;
        virtual T* end() = 0;
        virtual const T* end() const = 0;
        virtual const T* cbegin() const = 0;
        virtual const T* cend() const = 0;
        virtual bool empty() const = 0;
        virtual std::size_t size() const = 0;

        virtual bool resizable() const = 0;

        /// Overload iff resizable() returns true
        virtual void
        resize(bool)
        {
            throw std::runtime_error("Capsule cannot be resized");
        }
    };

    struct GenotypeCapsule : public Capsule<std::int8_t>
    {
        virtual ~GenotypeCapsule() = default;
        virtual std::size_t nsites() const = 0;
        virtual std::size_t nsam() const = 0;
        virtual std::size_t& nsites() = 0;
        virtual std::size_t& nsam() = 0;
        virtual std::size_t row_offset() const = 0;
        virtual std::size_t col_offset() const = 0;
        virtual std::size_t stride() const = 0;
        virtual std::unique_ptr<GenotypeCapsule> clone() const = 0;
        virtual std::int8_t& operator()(std::size_t, std::size_t) = 0;
        virtual const std::int8_t& operator()(std::size_t,
                                              std::size_t) const = 0;
    };

    struct PositionCapsule : public Capsule<double>
    {
        virtual ~PositionCapsule() = default;
        virtual std::size_t nsites() const = 0;
        virtual std::unique_ptr<PositionCapsule> clone() const = 0;
        virtual double& operator[](std::size_t) = 0;
        virtual const double& operator[](std::size_t) const = 0;
    };

} // namespace Sequence

#endif
