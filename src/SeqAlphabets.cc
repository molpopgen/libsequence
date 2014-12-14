//! \file src/SeqAlphabets.cc
#include <Sequence/SeqAlphabets.hpp>
#include <array>

namespace Sequence {
  const std::array<const char,17> dna_alphabet{ {'A','C','G','T',
	'R','Y','S','W',
	'K','M','B','D',
	'H','V','N','-','.'} };
}
