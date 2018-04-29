#ifndef SEQUENCE_VARIANT_MATRIX_INTERNAL_HPP__
#define SEQUENCE_VARIANT_MATRIX_INTERNAL_HPP__

#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <type_traits>
#include "col_view_iterator.hpp"

namespace Sequence
{
    namespace internal
    {
        template <typename T> struct row_view_
        /// Implementation details for Sequence::RowView and
        /// Sequence::ConstRowView
        /// \ingroup variantmatrix
        {
            static_assert(std::is_pointer<T>::value, "T must be pointer type");
            /// Pointer to row data
            T data;
            /// Data type
            using dtype = typename std::remove_pointer<T>::type;
            /// Number of elements in row.
            std::size_t row_size;

            row_view_(T data_, std::size_t row_size_)
                /// Constructor
                : data(data_), row_size(row_size_)
            {
            }
            inline dtype& operator[](const std::size_t i)
            /// Element access without range checking
            {
                return data[i];
            }
            inline const dtype& operator[](const std::size_t i) const
            /// Element access without range checking
            {
                return data[i];
            }
            inline dtype&
            at(const std::size_t i)
            /// Range-checked access
            {
                if (i >= row_size)
                    {
                        throw std::out_of_range("index out of range");
                    }
                return data[i];
            }
            inline const dtype&
            at(const std::size_t i) const
            /// Range-checked access
            {
                if (i >= row_size)
                    {
                        throw std::out_of_range("index out of range");
                    }
                return data[i];
            }
            std::size_t
            size() const
            /// Number of elements
            {
                return row_size;
            }

            using iterator = dtype*;
            using const_iterator = const dtype*;
            using reverse_iterator = std::reverse_iterator<iterator>;
            using const_reverse_iterator
                = std::reverse_iterator<const_iterator>;

            iterator
            begin()
            /// Get iterator to start of range
            {
                return data;
            }
            iterator
            end()
            /// Get iterator to end of range
            {
                return data + row_size;
            }
            const_iterator
            begin() const
            /// Get const iterator to start of range
            {
                return data;
            }
            const_iterator
            end() const
            /// Get const iterator to end of range
            {
                return data + row_size;
            }
            const_iterator
            cbegin() const
            /// Get const iterator to start of range
            {
                return this->begin();
            }
            const_iterator
            cend() const
            /// Get const iterator to end of range
            {
                return this->end();
            }

            // Reverse iterators
            reverse_iterator
            rbegin()
            /// Reverse iterator.  Points to start of reversed range.
            {
                return reverse_iterator(data + row_size);
            }
            reverse_iterator
            rend()
            /// Reverse iterator.  Points to end of reversed range.
            {
                return reverse_iterator(data);
            }
            const_reverse_iterator
            rbegin() const
            /// Const reverse iterator.  Points to start of reversed range.
            {
                return reverse_iterator(data + row_size);
            }
            const_reverse_iterator
            rend() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return reverse_iterator(data);
            }
            const_reverse_iterator
            crbegin() const
            /// Const reverse iterator.  Points to start of reversed range.
            {
                return this->rbegin();
            }
            const_reverse_iterator
            crend() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return this->rend();
            }

            std::vector<std::int8_t>
            copy() const
            /// Return copy of the view as std::vector<std::int8_t>
            {
                return std::vector<std::int8_t>(this->cbegin(), this->cend());
            }

            friend void
            swap(row_view_& a, row_view_& b)
            /// Allow swap via argument-dependent lookup, or "ADL"
            {
                if (a.size() != b.size())
                    {
                        throw std::invalid_argument(
                            "cannot swap row views of different size");
                    }
                auto bi = b.begin();
                for (auto ai = a.begin(); ai != a.end(); ++ai, ++bi)
                    {
                        std::swap(*ai, *bi);
                    }
            }
        };

        template <typename T> struct col_view_
        /// Implementation details for Sequence::ColView and
        /// Sequence::ConstColView
        /// \ingroup variantmatrix
        {
            static_assert(std::is_pointer<T>::value, "T must be pointer type");
            /// Pointer to column data
            T data;
            /// data type
            using dtype = typename std::remove_pointer<T>::type;

            /// data + col_end marks the end of the column data
            std::size_t col_end;
            /// Stride of the data in the column
            std::size_t stride;

            col_view_(T data_, std::size_t col_end_, std::size_t stride_)
                /// Constructor
                : data(data_), col_end(col_end_), stride(stride_)
            {
            }
            inline dtype& operator[](const std::size_t i)
            /// Element access without range checking.
            {
                return data[i * stride];
            }
            inline const dtype& operator[](const std::size_t i) const
            /// Element access without range checking.
            {
                return data[i * stride];
            }
            inline dtype&
            at(const std::size_t i)
            /// Range-checked access
            {
                if (i >= col_end / stride)
                    {
                        throw std::out_of_range("index out of range");
                    }
                return data[i * stride];
            }
            inline const dtype&
            at(const std::size_t i) const
            /// Range-checked access
            {
                if (i >= col_end / stride)
                    {
                        throw std::out_of_range("index out of range");
                    }
                return data[i * stride];
            }
            std::size_t
            size() const
            /// Number of elements
            {
                return col_end / stride;
            }

            using iterator = iterator_<dtype*>;
            using const_iterator = iterator_<const dtype*>;
            using reverse_iterator = std::reverse_iterator<iterator>;
            using const_reverse_iterator
                = std::reverse_iterator<const_iterator>;

            iterator
            begin()
            /// Get iterator to start of range
            {
                return iterator(
                    data,
                    static_cast<typename iterator::difference_type>(stride),
                    0);
            }
            iterator
            end()
            /// Get iterator to end of range
            {
                return iterator(
                    data,
                    static_cast<typename iterator::difference_type>(stride),
                    static_cast<typename iterator::difference_type>(col_end));
            }
            const_iterator
            begin() const
            /// Get const iterator to end of range
            {
                return const_iterator(
                    data,
                    static_cast<typename iterator::difference_type>(stride),
                    0);
            }
            const_iterator
            end() const
            /// Get const iterator to end of range
            {
                return const_iterator(
                    data,
                    static_cast<typename iterator::difference_type>(stride),
                    static_cast<typename iterator::difference_type>(col_end));
            }
            const_iterator
            cbegin() const
            /// Get const iterator to start of range
            {
                return this->begin();
            }
            const_iterator
            cend() const
            /// Get const iterator to end of range
            {
                return this->end();
            }

            // Reverse iterators
            reverse_iterator
            rbegin()
            /// Reverse iterator.  Points to start of reversed range.
            {
                return reverse_iterator(end());
            }
            reverse_iterator
            rend()
            /// Reverse iterator.  Points to end of reversed range.
            {
                return reverse_iterator(begin());
            }

            const_reverse_iterator
            rbegin() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return const_reverse_iterator(end());
            }
            const_reverse_iterator
            rend() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return const_reverse_iterator(begin());
            }
            const_reverse_iterator
            crbegin() const
            /// Const reverse iterator.  Points to start of reversed range.
            {
                return this->rbegin();
            }
            const_reverse_iterator
            crend() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return this->rend();
            }

            std::vector<std::int8_t>
            copy() const
            /// Return copy of the view as std::vector<std::int8_t>
            {
                return std::vector<std::int8_t>(this->cbegin(), this->cend());
            }

            friend void
            swap(col_view_& a, col_view_& b)
            /// Allow swap via argument-dependent lookup, or "ADL"
            {
                if (a.size() != b.size())
                    {
                        throw std::invalid_argument(
                            "cannot swap column views of different size");
                    }
                auto bi = b.begin();
                for (auto ai = a.begin(); ai != a.end(); ++ai, ++bi)
                    {
                        std::swap(*ai, *bi);
                    }
            }
        };
    }
}

#endif
