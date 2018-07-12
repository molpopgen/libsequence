#ifndef SEQUENCE_COL_VIEW_ITERATOR_HPP__
#define SEQUENCE_COL_VIEW_ITERATOR_HPP__

#include <iterator>
#include <stdexcept>

namespace Sequence
{
    template <typename POINTER> struct col_view_iterator
    /// \brief Iterator for column views.
    ///
    /// This is a C++11-compliant, random-access
    /// iterator
    /// \ingroup variantmatrix
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
        mutable POINTER start;
        /// Iterator data
        mutable POINTER data;

        /// Stride needed to increment/decrement
        difference_type stride;
        /// Offset w.r.to data
        difference_type offset;
        explicit col_view_iterator(POINTER data_, difference_type stride_,
                                   difference_type offset_)
            /// Constructor
            : start{ data_ }, data{ data_ }, stride{ stride_ }, offset{
                  offset_
              }
        {
        }

        col_view_iterator(const col_view_iterator&) = default;

        /// Get value pointed to
        reference operator*() { return *(data + offset); }

        /// Get value pointed to
        const reference operator*() const { return *(data + offset); }

        /// Access a reference n elements away
        reference operator[](difference_type n) { return *(*this + n); }
        /// Access a const reference n elements away
        const reference operator[](difference_type n) const
        {
            return *(*this + n);
        }

        col_view_iterator&
        operator=(const col_view_iterator& rhs)
        /// Assignment operator
        {
            this->start = rhs.start;
            this->data = rhs.data;
            this->stride = rhs.stride;
            this->offset = rhs.offset;
            return *this;
        }

        bool
        operator<=(const col_view_iterator rhs) const
        {
            return (this->data + this->offset) <= (rhs.data + rhs.offset);
        }
        bool
        operator<(const col_view_iterator rhs) const
        {
            return (this->data + this->offset) < (rhs.data + rhs.offset);
        }
        bool
        operator>(const col_view_iterator rhs) const
        {
            return !(*this <= rhs);
        }
        bool
        operator>=(const col_view_iterator rhs) const
        {
            return !(*this < rhs);
        }
        bool
        operator==(const col_view_iterator rhs) const
        {
            return !(*this < rhs) && !(*this > rhs);
        }
        bool
        operator!=(const col_view_iterator rhs) const
        {
            return !(*this == rhs);
        }
    };

    template <typename POINTER>
    col_view_iterator<POINTER> inline
    operator+(col_view_iterator<POINTER> i,
              typename col_view_iterator<POINTER>::difference_type d)
    /// Add to an iterator
    /// \ingroup variantmatrix
    {
        auto rv(i);
        rv.offset += d * rv.stride;
        return rv;
    }

    template <typename POINTER>
    col_view_iterator<POINTER> inline
    operator+(typename col_view_iterator<POINTER>::difference_type d,
              col_view_iterator<POINTER> i)
    /// Add to an iterator
    /// \ingroup variantmatrix
    {
        return i + d;
    }

    template <typename POINTER>
    inline col_view_iterator<POINTER>&
    operator++(col_view_iterator<POINTER>& i)
    /// Increment an iterator by one
    /// \ingroup variantmatrix
    {
        i.offset += i.stride;
        return i;
    }

    template <typename POINTER>
    inline col_view_iterator<POINTER>
    operator-(col_view_iterator<POINTER> i,
              typename col_view_iterator<POINTER>::difference_type d)
    /// Subtract from an iterator
    /// \ingroup variantmatrix
    {
        auto rv(i);
        rv.offset -= d * rv.stride;
        return rv;
    }

    template <typename POINTER>
    col_view_iterator<POINTER> inline
    operator-(typename col_view_iterator<POINTER>::difference_type d,
              col_view_iterator<POINTER> i)
    /// Subtract from an iterator
    /// \ingroup variantmatrix
    {
        return i - d;
    }

    template <typename POINTER>
    inline col_view_iterator<POINTER>&
    operator--(col_view_iterator<POINTER>& i)
    /// Decrement iterator by one
    /// \ingroup variantmatrix
    {
        i.offset -= i.stride;
        return i;
    }

    template <typename POINTER>
    inline col_view_iterator<POINTER>&
    operator+=(col_view_iterator<POINTER>& i,
               typename col_view_iterator<POINTER>::difference_type d)
    /// Increment an iterator
    /// \ingroup variantmatrix
    {
        i.offset += d * i.stride;
        return i;
    }

    template <typename POINTER>
    inline col_view_iterator<POINTER>&
    operator-=(col_view_iterator<POINTER>& i,
               typename col_view_iterator<POINTER>::difference_type d)
    /// Decrement an iterator
    /// \ingroup variantmatrix
    {
        i.offset -= d * i.stride;
        return i;
    }

    template <typename POINTER>
    inline typename col_view_iterator<POINTER>::difference_type
    operator-(col_view_iterator<POINTER> i, col_view_iterator<POINTER> j)
    /// Distance between two iterators
    /// \ingroup variantmatrix
    {
        if (i.start != j.start)
            {
                throw std::invalid_argument("attempt to subtract "
                                            "iterators from different "
                                            "columns");
            }
        return (i.offset - j.offset) / i.stride;
    }
}

#endif
