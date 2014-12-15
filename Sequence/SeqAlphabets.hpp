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

  /*
    @brief test if character is part of Sequence::dna_alphabet
    @param ch Character to test
    \return true if ch is in Sequence::dna_alphabet, false otherwise
    \note case-insensitive via std::toupper
  */
  bool isDNA( const char & ch);
}

#endif
