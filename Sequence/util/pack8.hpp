#ifndef __SEQUENCE_UTIL_PACK8_HPP__
#define __SEQUENCE_UTIL_PACK8_HPP__
#include <Sequence/SeqAlphabets.hpp>
#include <cctype>
#include <string>
#include <vector>
#include <array>

namespace Sequence
{
  namespace pack8
  {
    using itype = std::uint8_t;
    using vtype = std::vector<itype>;

    vtype dna2vtype( const std::string & s, const Sequence::alphabet_t & a );
    vtype dna2vtype( std::string & s, const Sequence::alphabet_t & a );
    std::string vtype2dna( const vtype & v, const Sequence::alphabet_t & a, const std::string::size_type & len );
    std::string vtype2dna( vtype & v, const Sequence::alphabet_t & a, const std::string::size_type & len );
  }
}
#endif
