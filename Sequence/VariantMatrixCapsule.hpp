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
        virtual T& operator[](std::size_t) = 0;
        virtual const T& operator[](std::size_t) const = 0;
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
    };

    struct GenotypeCapsule : public Capsule<std::int8_t>
    {
        virtual ~GenotypeCapsule() = default;
        virtual std::size_t nsites() const = 0;
        virtual std::size_t nsam() const = 0;
        virtual std::size_t& nsites() = 0;
        virtual std::size_t& nsam() = 0;
        virtual std::unique_ptr<GenotypeCapsule> clone() const = 0;
    };

    struct PositionCapsule : public Capsule<double>
    {
        virtual ~PositionCapsule() = default;
        virtual std::size_t nsites() const = 0;
        virtual std::unique_ptr<PositionCapsule> clone() const = 0;
    };

    class VectorGenotypeCapsule : public GenotypeCapsule
    {
      private:
        std::vector<std::int8_t> buffer;
        std::size_t nsites_, nsam_;

        std::size_t
        fill_nsam(std::size_t num_rows)
        {
            return ((num_rows > 0) ? buffer.size() / num_rows : 0);
        }

      public:
        template <typename T>
        VectorGenotypeCapsule(T&& t, std::size_t num_rows)
            : buffer(std::forward<T>(t)), nsites_(num_rows),
              nsam_(fill_nsam(num_rows))
        {
        }

        std::size_t&
        nsites()
        {
            return nsites_;
        }

        std::size_t&
        nsam()
        {
            return nsam_;
        }

        std::size_t
        nsites() const
        {
            return nsites_;
        }

        std::size_t
        nsam() const
        {
            return nsam_;
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

        const std::int8_t*
        cdata() const final
        {
            return buffer.data();
        }

        std::unique_ptr<GenotypeCapsule>
        clone() const final
        {
            return std::unique_ptr<GenotypeCapsule>(
                new VectorGenotypeCapsule(this->buffer, this->nsites_));
        }

        std::int8_t*
        begin() final
        {
            return buffer.data();
        }

        const std::int8_t*
        begin() const final
        {
            return buffer.data();
        }

        std::int8_t*
        end() final
        {
            return buffer.data() + buffer.size();
        }

        const std::int8_t*
        end() const final
        {
            return buffer.data() + buffer.size();
        }

        const std::int8_t*
        cbegin() const final
        {
            return begin();
        }

        const std::int8_t*
        cend() const final
        {
            return end();
        }

        bool
        empty() const final
        {
            return buffer.empty();
        }

        std::size_t
        size() const final
        {
            return buffer.size();
        }
    };

    class VectorPositionCapsule : public PositionCapsule
    {
      private:
        std::vector<double> buffer;

      public:
        template <typename T>
        explicit VectorPositionCapsule(T&& t) : buffer(std::forward<T>(t))
        {
        }
        double& operator[](std::size_t i) { return buffer[i]; }

        const double& operator[](std::size_t i) const { return buffer[i]; }

        double*
        data() final
        {
            return buffer.data();
        }

        const double*
        data() const final
        {
            return buffer.data();
        }

        const double*
        cdata() const final
        {
            return buffer.data();
        }

        std::unique_ptr<PositionCapsule>
        clone() const final
        {
            return std::unique_ptr<PositionCapsule>(
                new VectorPositionCapsule(*this));
        }

        double*
        begin() final
        {
            return buffer.data();
        }

        const double*
        begin() const final
        {
            return buffer.data();
        }

        double*
        end() final
        {
            return buffer.data() + buffer.size();
        }

        const double*
        end() const final
        {
            return buffer.data() + buffer.size();
        }

        const double*
        cbegin() const final
        {
            return begin();
        }

        const double*
        cend() const final
        {
            return end();
        }

        bool
        empty() const final
        {
            return buffer.empty();
        }

        std::size_t
        size() const final
        {
            return buffer.size();
        }

        std::size_t
        nsites() const
        {
            return buffer.size();
        }
    };
} // namespace Sequence

#endif
