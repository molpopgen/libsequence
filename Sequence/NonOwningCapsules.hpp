#ifndef NONOWNING_CAPSULES_HPP
#define NONOWNING_CAPSULES_HPP

#include "VariantMatrixCapsule.hpp"
#include <vector>

namespace Sequence
{
    class NonOwningGenotypeCapsule : public GenotypeCapsule
    {
      private:
        struct nodelete
        {
            void
            operator()(const std::int8_t*)
            {
            }
        };
        std::unique_ptr<const std::int8_t[], nodelete> buffer;
        std::size_t nsites_, nsam_, k1, k2, tda;

      public:
        NonOwningGenotypeCapsule(const std::int8_t* data, std::size_t nrow,
                                 std::size_t ncol, std::size_t trailing)
            : buffer(data, nodelete()), nsites_(nrow), nsam_(ncol), k1(0),
              k2(0), tda(trailing)
        {
        }

        NonOwningGenotypeCapsule(const std::int8_t* data, std::size_t nrow,
                                 std::size_t ncol, std::size_t row_offset,
                                 std::size_t column_offset, std::size_t trailing)
            : buffer(data, nodelete()), nsites_(nrow), nsam_(ncol),
              k1(row_offset), k2(column_offset), tda(trailing)
        {
        }

        std::size_t& nsites();

        std::size_t& nsam();

        std::size_t nsites() const;

        std::size_t nsam() const;

        [[noreturn]] std::int8_t& operator()(std::size_t, std::size_t);

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
    };

    class NonOwningPositionCapsule : public PositionCapsule
    {
      private:
        struct nodelete
        {
            void
            operator()(const double*)
            {
            }
        };
        std::unique_ptr<const double[], nodelete> buffer;
        std::size_t current_size;

      public:
        NonOwningPositionCapsule(const double* data, const std::size_t nsites)
            : buffer(data, nodelete()), current_size(nsites)
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
    };
} // namespace Sequence
#endif

