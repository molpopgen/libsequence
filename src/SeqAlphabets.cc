//! \file src/SeqAlphabets.cc
#include <Sequence/SeqAlphabets.hpp>
#include <algorithm>
#include <cctype>
#include <array>

namespace Sequence {
  const std::array<const char,17> dna_alphabet{ {'A','C','G','T',
	'R','Y','S','W',
	'K','M','B','D',
	'H','V','N','-','.'} };

  bool isDNA( const char & ch) 
  {
    return std::find( dna_alphabet.begin(),
		      dna_alphabet.end(),
		      std::toupper(ch) ) != dna_alphabet.end();
  }
}
