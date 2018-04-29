#ifndef SEQUENCE_COL_VIEW_ITERATOR_HPP__
#define SEQUENCE_COL_VIEW_ITERATOR_HPP__

#include <iterator>

namespace Sequence
{
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
        using value_type = typename std::iterator_traits<POINTER>::value_type;
        /// Reference type
        using reference = typename std::iterator_traits<POINTER>::reference;
        /// Pointer type
        using pointer = POINTER;
        /// Iterator category
        using iterator_category =
            typename std::iterator_traits<POINTER>::iterator_category;

        /// The start of a view
        /// Used to ensure that two
        /// column iterators refer to the
        /// same column
        const POINTER start;
        /// Iterator data
        mutable POINTER data;

        /// Stride needed to increment/decrement
        difference_type stride;
        /// Offset w.r.to data
        difference_type offset;
        explicit iterator_(POINTER data_, difference_type stride_,
                           difference_type offset_)
            /// Constructor
            : start{ data_ }, data{ data_ }, stride{ stride_ }, offset{
                  offset_
              }
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

        difference_type
        operator-(iterator_ i)
        {
            if (this->start != i.start)
                {
                    throw std::invalid_argument("attempt to subtract "
                                                "iterators from different "
                                                "columns");
                }
            return (this->offset - i.offset) / this->stride;
        }

        bool
        operator<=(const iterator_ rhs) const
        {
            return (this->data + this->offset) <= (rhs.data + rhs.offset);
        }
        bool
        operator<(const iterator_ rhs) const
        {
            return (this->data + this->offset) < (rhs.data + rhs.offset);
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
}

#endif
