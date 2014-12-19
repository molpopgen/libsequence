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
    using itype = std::int8_t;
    //! vectors of itype
    using vtype = std::vector<itype>;

    vtype dna2vtype( const std::string & s );
    vtype dna2vtype( std::string & s );
    vtype dna2vtype( const Sequence::Seq & s );
    vtype dna2vtype( Sequence::Seq & s );
    std::string vtype2dna( const vtype & v );
    std::string vtype2dna( vtype & v );
  }
}

#endif
