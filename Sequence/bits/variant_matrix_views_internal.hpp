#ifndef SEQUENCE_VARIANT_MATRIX_INTERNAL_HPP__
#define SEQUENCE_VARIANT_MATRIX_INTERNAL_HPP__

#include <iterator>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <type_traits>

namespace Sequence
{
    namespace internal
    {
        template <typename T> struct row_view_
        /// Implementation details for Sequence::RowView and
        /// Sequence::ConstRowView
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

            template <typename POINTER> struct iterator_
            /// Iterator for column views.
            /// This is a C++11 standard compliant iterator.
            {
                static_assert(std::is_pointer<POINTER>::value,
                              "iterator must wrap a pointer type");
                /// Difference type
                using difference_type =
                    typename std::iterator_traits<POINTER>::difference_type;
                /// Value type
                using value_type =
                    typename std::iterator_traits<POINTER>::value_type;
                /// Reference type
                using reference =
                    typename std::iterator_traits<POINTER>::reference;
                /// Pointer type
                using pointer = POINTER;
                /// Iterator category
                using iterator_category =
                    typename std::iterator_traits<POINTER>::iterator_category;

                /// Iterator data
                mutable POINTER data;

                /// Stride needed to increment/decrement
                difference_type stride;
                /// Offset w.r.to data
                difference_type offset;
                explicit iterator_(POINTER data_, difference_type stride_,
                                   difference_type offset_)
                    : data{ data_ }, stride{ stride_ }, offset{ offset_ }
                /// Constructor
                {
                }

                /// Get value pointed to
                reference operator*() { return *(data + offset); }

                /// Get value pointed to
                const reference operator*() const { return *(data + offset); }

                iterator_&
                operator=(const iterator_& rhs)
                /// Assignment operator
                {
                    this->data = rhs.data;
                    this->stride = rhs.stride;
                    this->offset = rhs.offset;
                }
                iterator_&
                operator+(difference_type d)
                {
                    offset += d * stride;
                    return *this;
                }
                iterator_& operator++() { return *this + 1; }
                iterator_&
                operator+=(difference_type d)
                {
                    return *this + d;
                }

                iterator_&
                operator-(difference_type d)
                {
                    return *this + -d;
                }
                iterator_& operator--() { return *this - 1; }
                iterator_&
                operator-=(difference_type d)
                {
                    return *this - d;
                }
                bool
                operator<=(const iterator_ rhs) const
                {
                    return (this->data + this->offset)
                           <= (rhs.data + rhs.offset);
                }
                bool
                operator<(const iterator_ rhs) const
                {
                    return (this->data + this->offset)
                           < (rhs.data + rhs.offset);
                }
                bool
                operator>(const iterator_ rhs) const
                {
                    return !(*this <= rhs);
                }
                bool
                operator>=(const iterator_ rhs) const
                {
                    return !(*this < rhs);
                }
                bool
                operator==(const iterator_ rhs) const
                {
                    return !(*this < rhs) && !(*this > rhs);
                }
                bool
                operator!=(const iterator_ rhs) const
                {
                    return !(*this == rhs);
                }
            };

            using iterator = iterator_<dtype*>;
            using const_iterator = iterator_<const dtype*>;
            using reverse_iterator = std::reverse_iterator<iterator>;
            using const_reverse_iterator
                = std::reverse_iterator<const_iterator>;

            iterator
            begin()
            /// Get iterator to start of range
            {
                return iterator(data, stride, 0);
            }
            iterator
            end()
            /// Get iterator to end of range
            {
                return iterator(data, stride, col_end);
            }
            const_iterator
            begin() const
            /// Get const iterator to end of range
            {
                return const_iterator(data, stride, 0);
            }
            const_iterator
            end() const
            /// Get const iterator to end of range
            {
                return const_iterator(data, stride, col_end);
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
                return reverse_iterator(iterator(data, stride, col_end));
            }
            reverse_iterator
            rend()
            /// Reverse iterator.  Points to end of reversed range.
            {
                return reverse_iterator(iterator(data, stride, 0));
            }

            const_reverse_iterator
            rbegin() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return const_reverse_iterator(
                    const_iterator(data, stride, col_end));
            }
            const_reverse_iterator
            rend() const
            /// Const reverse iterator.  Points to end of reversed range.
            {
                return const_reverse_iterator(const_terator(data, stride, 0));
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
