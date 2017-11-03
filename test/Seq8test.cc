/*! \file Seq8test.cc @brief Unit tests for Sequence/Seq.hpp */
#define BOOST_TEST_MODULE Seq8test

#include <Sequence/Seq8.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iostream>
BOOST_AUTO_TEST_CASE( create_1 )
{
  std::string seq("ATGCN");
  Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet);

  BOOST_CHECK_EQUAL(seq.size(),s8.size());

  for( unsigned i = 0 ; i < seq.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( seq[i],s8[i] );
    }

  BOOST_CHECK_EQUAL( seq, s8.unpack() );
}

BOOST_AUTO_TEST_CASE( create_2 )
{
  std::string seq("ATGCNX");
  BOOST_CHECK_THROW(Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet),
		    Sequence::SeqException);

}

BOOST_AUTO_TEST_CASE( create_3 )
{
  std::string seq("ATGCNW"); //valid DNA, invalid POLY chars
  BOOST_CHECK_THROW(Sequence::Seq8 s8(std::cref(seq), Sequence::dna_poly_alphabet),
		    Sequence::SeqException);

}

BOOST_AUTO_TEST_CASE( assign_1 )
{
  std::string seq("ATGCN");
  Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet),
    s8_2(s8);

  BOOST_CHECK(s8==s8_2);
}

BOOST_AUTO_TEST_CASE( assign_2 )
{
  std::string seq("ATGCN");
  Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet),
    s8_2 = s8;

  BOOST_CHECK(s8==s8_2);
}

BOOST_AUTO_TEST_CASE( assign_3 )
{
  std::string seq("ATGCN");
  Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet),
    s8_2(std::move(s8));

  BOOST_CHECK_EQUAL( s8.second.size(), 0 );

  BOOST_CHECK_EQUAL(seq.size(),s8_2.size());

  for( unsigned i = 0 ; i < seq.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( seq[i],s8_2[i] );
    }
}

BOOST_AUTO_TEST_CASE( assign_4 )
{
  std::string seq("ATGCN");
  Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet),
    s8_2 = std::move(s8);

  BOOST_CHECK_EQUAL( s8.second.size(), 0 );

  BOOST_CHECK_EQUAL(seq.size(),s8_2.size());

  for( unsigned i = 0 ; i < seq.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( seq[i],s8_2[i] );
    }
}


BOOST_AUTO_TEST_CASE( IO_1 )
{
  std::string seq("ATGCN");
  Sequence::Seq8 s8(std::cref(seq), Sequence::dna_alphabet),
    s8_2;

  std::ostringstream o;
  o << s8;
  std::istringstream i(o.str());
  i >> s8_2;
  BOOST_CHECK( s8.first == s8_2.first );
  BOOST_CHECK( s8.second == s8_2.second );
  BOOST_CHECK_EQUAL( s8.unpack() , s8_2.unpack() );
}

