/*! \file alphabets.cc @brief Unit tests for Sequence/SeqAlphabets.hpp */

#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/Fasta.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>
BOOST_AUTO_TEST_SUITE(AlphabetTest)

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

//Test of dna_poly_alphabet
BOOST_AUTO_TEST_CASE( dna_poly_alphabet_1 )
{
  BOOST_REQUIRE( std::find( Sequence::dna_poly_alphabet.begin(),
			    Sequence::dna_poly_alphabet.end(),'\0' ) != Sequence::dna_poly_alphabet.end() );
}

BOOST_AUTO_TEST_CASE( dna_poly_alphabet_2 )
{
  BOOST_CHECK_EQUAL( Sequence::POLYEOS, 8 );
}

BOOST_AUTO_TEST_CASE( dna_poly_alphabet_3 )
{
  for ( auto c : {'A','C','G','T','N','0','1','-'} )
    {
      BOOST_CHECK( std::distance( Sequence::dna_poly_alphabet.begin(),
				  std::find(Sequence::dna_poly_alphabet.begin(),
					    Sequence::dna_poly_alphabet.end(),c) ) < Sequence::POLYEOS );
    }
}

BOOST_AUTO_TEST_CASE( dna_poly_alphabet_4 )
{
  for ( auto c : {'a','c','g','t','n','W','K'} )
    {
      BOOST_CHECK( std::distance( Sequence::dna_poly_alphabet.begin(),
				  std::find(Sequence::dna_poly_alphabet.begin(),
					    Sequence::dna_poly_alphabet.end(),c) ) >= Sequence::POLYEOS );
    }
}

BOOST_AUTO_TEST_CASE( dna_poly_alphabet_5 )
{
  for ( auto c : {'a','c','g','t','n','W','K'} )
    {
      BOOST_CHECK( std::distance( Sequence::dna_poly_alphabet.begin(),
				  std::find(Sequence::dna_poly_alphabet.begin(),
					    Sequence::dna_poly_alphabet.end(),c) ) >= Sequence::NOTPOLYCHAR );
    }
}
BOOST_AUTO_TEST_SUITE_END()
