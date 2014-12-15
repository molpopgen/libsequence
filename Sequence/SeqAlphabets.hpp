//! \file Sequence/SeqAlphabets.hpp
#ifndef __SEQUENCE_SEQALPHABETS_HPP__
#define __SEQUENCE_SEQALPHABETS_HPP__

#include <array>

namespace Sequence {
  /*!
    @brief Alphabet for DNA sequences
    Valid DNA characters.  Upper-case only.  
    Both . and - are accepted as gap characters
    \note http://www.bioinformatics.org/sms/iupac.html, excluding U
  */
  extern const std::array<const char,17> dna_alphabet;
}

#endif
