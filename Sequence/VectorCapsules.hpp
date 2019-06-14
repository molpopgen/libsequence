#ifndef VECTOR_CAPSULES_HPP
#define VECTOR_CAPSULES_HPP

#include "VariantMatrixCapsule.hpp"
#include <vector>

namespace Sequence
{
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

        std::size_t& nsites();

        std::size_t& nsam();

        std::size_t nsites() const;

        std::size_t nsam() const;

        std::int8_t& operator()(std::size_t, std::size_t);

        const std::int8_t& operator()(std::size_t, std::size_t) const;

        std::int8_t* data() final;

        const std::int8_t* data() const final;

        const std::int8_t* cdata() const final;

        std::unique_ptr<GenotypeCapsule> clone() const final;

        std::int8_t* begin() final;

        const std::int8_t* begin() const final;

        std::int8_t* end() final;

        const std::int8_t* end() const final;

        const std::int8_t* cbegin() const final;

        const std::int8_t* cend() const final;

        bool empty() const final;

        std::size_t size() const final;

        bool resizable() const final;

        void resize(bool) final;
    };

    class VectorPositionCapsule : public PositionCapsule
    {
      private:
        std::vector<double> buffer;
        std::size_t current_size;

      public:
        template <typename T>
        explicit VectorPositionCapsule(T&& t)
            : buffer(std::forward<T>(t)), current_size(buffer.size())
        {
        }

        double& operator[](std::size_t);

        const double& operator[](std::size_t) const;

        double* data() final;

        const double* data() const final;

        const double* cdata() const final;

        std::unique_ptr<PositionCapsule> clone() const final;

        double* begin() final;

        const double* begin() const final;

        double* end() final;

        const double* end() const final;

        const double* cbegin() const final;

        const double* cend() const final;

        bool empty() const final;

        std::size_t size() const final;

        std::size_t nsites() const;

        bool resizable() const final;

        void resize(bool) final;
    };
} // namespace Sequence
#endif
