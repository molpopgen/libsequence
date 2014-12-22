//! \file Sequence/util/pack8.hpp @brief namespace Sequence::pack8
#ifndef __SEQUENCE_UTIL_PACK8_HPP__
#define __SEQUENCE_UTIL_PACK8_HPP__
#include <Sequence/SeqAlphabets.hpp>
#include <cctype>
#include <string>
#include <vector>
#include <array>

/*! \namespace Sequence::pack8
  @brief routines for performing 4-bit encoding of DNA sequences into 8-bit integers
 */
namespace Sequence
{
  namespace pack8
  {
    //! The integer type
    using itype = std::uint8_t;
    //! The sequence type
    using vtype = std::vector<itype>;

    /*!
      Pack a sequence using an alphabet
    */
    vtype dna2vtype( const std::string & s, const Sequence::alphabet_t & a );
    /*!
      Pack a sequence using an alphabet
      \note s will be cleared
    */
    vtype dna2vtype( std::string & s, const Sequence::alphabet_t & a );
    /*!
      Unpack a sequence using an alphabet
    */
    std::string vtype2dna( const vtype & v, const Sequence::alphabet_t & a, const std::string::size_type & len );
    /*!
      Unpack a sequence using an alphabet
      \note v will be cleared
    */
    std::string vtype2dna( vtype & v, const Sequence::alphabet_t & a, const std::string::size_type & len );
  }
}
#endif
