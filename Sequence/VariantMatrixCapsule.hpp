#ifndef SEQUENCE_VARIANT_MATRIX_CAPSULE
#define SEQUENCE_VARIANT_MATRIX_CAPSULE

#include <cstdint>
#include <memory>

namespace Sequence
{
    struct VariantMatrixCapsule
    {
        virtual ~VariantMatrixCapsule() = default;
        virtual std::int8_t& operator[](std::size_t) = 0;
        virtual const std::int8_t& operator[](std::size_t) const = 0;
        // Following two may not be needed
        //virtual std::int8_t& get(std::size_t site, std::size_t sample) = 0;
        //virtual const std::int8_t& get(std::size_t site,
        //                               std::size_t sample) const = 0;
        virtual std::int8_t* data() = 0;
        virtual const std::int8_t* data() const = 0;
        virtual std::unique_ptr<VariantMatrixCapsule> clone() const = 0;
    };

    class VectorCapsule : public VariantMatrixCapsule
    {
      private:
        std::vector<std::int8_t> buffer;

      public:
        template <typename T>
        explicit VectorCapsule(T&& t) : buffer(std::forward<T>(t))
        {
        }

        std::int8_t& operator[](std::size_t i) { return buffer[i]; }

        const std::int8_t& operator[](std::size_t i) const
        {
            return buffer[i];
        }

        std::int8_t*
        data() final
        {
            return buffer.data();
        }

        const std::int8_t*
        data() const final
        {
            return buffer.data();
        }

        std::unique_ptr<VariantMatrixCapsule>
        clone() const final
        {
            return std::unique_ptr<VariantMatrixCapsule>(
                new VectorCapsule(this->buffer));
        }
    };
} // namespace Sequence

#endif
