//! \file Sequence/SeqAlphabets.hpp
#ifndef __SEQUENCE_SEQALPHABETS_HPP__
#define __SEQUENCE_SEQALPHABETS_HPP__

#include <array>
#include <functional>

namespace Sequence {
  /*!
    @brief Alphabet for DNA sequences
    Valid DNA characters.  Upper-case only.  
    Both . and - are accepted as gap characters
    \note http://www.bioinformatics.org/sms/iupac.html, excluding U
  */
  extern const std::array<const char,17> dna_alphabet;

  /*!
    @brief test if character is part of Sequence::dna_alphabet
    @param ch Character to test
    \return true if ch is in Sequence::dna_alphabet, false otherwise
    \note case-insensitive via std::toupper
  */
  bool isDNA( const char & ch);

  /*!
    \struct ambiguousNucleotide Sequence/SeqAlphabets.hpp
    @brief Tests if a character is in the set A,G,C,T
  */
  struct ambiguousNucleotide : public std::unary_function<char,bool>
  {
    /*!
      \return true if c is not A,G,C, or T, false otherwise
      \note Case-insensitive
    */
    bool operator()(const char & c) const;
  };

  /*!
    \struct invalidPolyChar Sequence/SeqAlphabets.hpp
    @brief This functor can be used to determine
    if a range contains characters that
    the SNP analysis routines in this
    library cannot handle gracefully
  */
  struct invalidPolyChar : public std::unary_function<char,bool>
  {
    /*!
      \return true if c is not in the set {A,G,C,T,N,-,.,0,1}, false otherwise. The period
      (.) can be used as an "identical to the 1st seq in a file" character, so 
      should be considered valid
      \note Case-insensitive
    */
    bool operator()(const char & nucleotide) const;
  };
}

#endif
