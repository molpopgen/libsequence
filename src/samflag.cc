#include <Sequence/samflag.hpp>

namespace Sequence
{
  samflag::operator int()
  /*!
    typecast to integer
    \return The value of the SAM flag
  */
  {
    return flag;
  }

  samflag::operator int() const
  /*!
    typecast to integer
    \return The value of the SAM flag
  */
  {
    return flag;
  }

  std::istream & samflag::read( std::istream & i)
  /*!
    Reads from an istream.  Generally only
    called via operator>>
  */
  {
    i >> flag;
    process_bits();
    return i;
  }

  void samflag::process_bits()
  /*!
    Takes the integer value of a flag
    and assigns to the boolean members of this class
  */
  {
    is_paired = (flag & sambits::is_paired);
    is_proper_pair = (flag & sambits::is_proper_pair);
    query_unmapped = (flag & sambits::query_unmapped);
    mate_unmapped = (flag & sambits::mate_unmapped);
    qstrand = (flag & sambits::qstrand);
    mstrand = (flag & sambits::mstrand);
    first_read = (flag & sambits::first_read);
    second_read = (flag & sambits::second_read);
    not_primary = (flag & sambits::not_primary);
    qcfail = (flag & sambits::qcfail);
    duplicate = (flag & sambits::duplicate);
  }

  std::ostream & operator<<(std::ostream & o, const samflag & s)
  /*!
    \ingroup operators
    Write a Sequence::samflag
  */
  {
    o << s.flag;
    return o;
  }

  std::istream & operator>>(std::istream & i, samflag & s)
  /*!
    \ingroup operators
    Read a Sequence::samflag
  */
  {
    return s.read(i);
  }

  //variable documentation
  /*!
    \var Sequence::samflag::duplicate
    True if (flag & Sequence::sambits::duplicate)
  */
  /*!
    \var Sequence::samflag::qcfail
    True if (flag & Sequence::sambits::qcfail)
  */
  /*!
    \var Sequence::samflag::not_primary
    True if (flag & Sequence::sambits::not_primary)
  */
  /*!
    \var Sequence::samflag::second_read
    True if (flag & Sequence::sambits::second_read)
  */
  /*!
    \var Sequence::samflag::first_read
    True if (flag & Sequence::sambits::first_read)
  */
  /*!
    \var Sequence::samflag::mstrand
    True if (flag & Sequence::sambits::mstrand)
  */
  /*!
    \var Sequence::samflag::qstrand
    True if (flag & Sequence::sambits::qstrand)
  */
  /*!
    \var Sequence::samflag::is_paired
    True if (flag & Sequence::sambits::is_paired)
  */
  /*!
    \var Sequence::samflag::is_proper_pair
    True if (flag & Sequence::sambits::is_proper_pair)
  */
  /*!
    \var Sequence::samflag::query_unmapped
    True if (flag & Sequence::sambits::query_unmapped)
  */
  /*!
    \var Sequence::samflag::mate_unmapped
    True if (flag & Sequence::sambits::mate_unmapped)
  */
}
