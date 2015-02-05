//! \file src/SeqAlphabets.cc
#include <Sequence/SeqAlphabets.hpp>
#include <algorithm>
#include <cctype>
#include <array>

namespace Sequence {
  const alphabet_t dna_alphabet{ {'A','C','G','T',
	'R','Y','S','W',
	'K','M','B','D',
	'H','V','N','-'} };

  const alphabet_t dna_poly_alphabet{ {'A','C','G','T', //0-3
	'0','1','-','N', //4-7
	'\0', //8
	} };

  const alphabet_t::size_type NOTPOLYCHAR = dna_poly_alphabet.size();
  
  const alphabet_t::size_type POLYEOS = alphabet_t::size_type( std::distance(dna_poly_alphabet.begin(),
									     std::find(dna_poly_alphabet.begin(),
										       dna_poly_alphabet.end(),
										       '\0')
									     ) );
  bool isDNA( const char & ch) 
  {
    return std::find( dna_alphabet.begin(),
		      dna_alphabet.end(),
		      std::toupper(ch) ) != dna_alphabet.end();
  }

  bool ambiguousNucleotide::operator()(const char & c) const
  {
    return std::distance( dna_alphabet.begin(),
			  std::find(dna_alphabet.begin(),
				    dna_alphabet.end(),
				    std::toupper(c)) ) > 3;
    /*
    const char ch = char(std::toupper(c));
    return (ch != 'A' &&
	    ch != 'G' &&
	    ch != 'T' &&
	    ch != 'C' );
    */
  }
  
  bool invalidPolyChar::operator()(const char & nucleotide) const
    {
      auto itr = std::find(dna_alphabet.begin(),
			   dna_alphabet.end(),
			   std::toupper(nucleotide));
      if(itr == dna_alphabet.end()) return 1;
      auto d = std::distance( dna_alphabet.begin(),
			      itr );
      return ( d > 3 && d < 14 ); 
    }
}
