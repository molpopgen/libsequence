/*! \file alphabets.cc @brief Unit tests for Sequence/SeqAlphabets.hpp */
#define BOOST_TEST_MODULE alphabets
#define BOOST_TEST_DYN_LINK 

#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/Fasta.hpp>
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

BOOST_AUTO_TEST_CASE( check_isDNA_1 )
{
  for (auto c : Sequence::dna_alphabet )
    {
      BOOST_REQUIRE( Sequence::isDNA(c) );
    }
}

BOOST_AUTO_TEST_CASE( check_isDNA_2 )
{
  Sequence::Fasta f = { "name","ATGCZAGC" };  //Z is a non-DNA character
  auto itr = std::find_if( f.begin(),f.end(),
			   [](const char & __ch) {
			     return !Sequence::isDNA(__ch);
			   } );
  BOOST_REQUIRE_EQUAL( std::distance(f.begin(),itr),4 );

  f.second.erase( std::remove_if(f.begin(),
				 f.end(),
				 [](const char & __ch) {
				   return !Sequence::isDNA(__ch);
				 }), f.second.end() );
  BOOST_REQUIRE_EQUAL(f.second,"ATGCAGC");
}
