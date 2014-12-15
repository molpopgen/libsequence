/*! \file alphabets.cc @brief Unit tests for Sequence/SeqAlphabets.hpp */
#define BOOST_TEST_MODULE alphabets
#define BOOST_TEST_DYN_LINK 

#include <Sequence/SeqAlphabets.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>

BOOST_AUTO_TEST_CASE( check_dna_alphabet )
{
  for ( auto c : {'A','G','C','T'} )
    {
      BOOST_REQUIRE( std::distance(Sequence::dna_alphabet.begin(),
				   std::find( Sequence::dna_alphabet.begin(),
					      Sequence::dna_alphabet.end(), c ) ) < 4 );
    }
}

BOOST_AUTO_TEST_CASE( check_isDNA )
{
  for (auto c : Sequence::dna_alphabet )
    {
      BOOST_REQUIRE( Sequence::isDNA(c) );
    }
}
