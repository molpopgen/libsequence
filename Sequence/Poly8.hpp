#ifndef __SEQUENCE_POlY8_HPP__
#define __SEQUENCE_POlY8_HPP__

#include <vector>
#include <cstdint>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/Seq.hpp>

namespace Sequence
{
  /*!
    \namespace Sequence::poly8
    @brief Converting std::strings to std::vector<std::int8_t> using Sequence::dna_poly_alphabet
  */
  namespace poly8
  {
    //! Integer type
    using itype = std::uint8_t;
    //! vectors of itype
    using vtype = std::vector<itype>;

    /*!
      \return A vtype made from s
     */
    vtype dna2vtype( const std::string & s );
    /*!
      \return A vtype made from s
      \note Contents of s are cleared
     */
    vtype dna2vtype( std::string & s );

    /*!
      \return A vtype made from s
      \note Contents of s are cleared
     */
    vtype dna2vtype( std::string && s );

    /*!
      \return A vtype made from s
     */
    vtype dna2vtype( const Sequence::Seq & s );

    /*!
      \return A vtype made from s
      \note Contents of s are cleared
     */
    vtype dna2vtype( Sequence::Seq & s );

    /*!
      \return A vtype made from s
      \note Contents of s are cleared
     */
    vtype dna2vtype( Sequence::Seq && s );

    /*!
      \return A std::string made from v
    */
    std::string vtype2dna( const vtype & v );

    /*!
      \return A std::string made from v
      \note Contents of v are cleared
    */
    std::string vtype2dna( vtype & v );

    /*!
      \return A std::string made from v
      \note Contents of v are cleared
    */
    std::string vtype2dna( vtype && v );
  }
}

#endif
