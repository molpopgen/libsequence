/*! \file pack8test.cc @brief Unit tests for Sequence/util/pack8.hpp and namespace Sequence::pack8*/
#define BOOST_TEST_MODULE pack8test
#define BOOST_TEST_DYN_LINK 

#include <Sequence/util/pack8.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>

BOOST_AUTO_TEST_CASE( even_dna_string )
{
  std::string seq("ATGCN");
  Sequence::pack8::vtype v = Sequence::pack8::dna2vtype(std::cref(seq), Sequence::dna_alphabet);
  std::string seq2 = Sequence::pack8::vtype2dna( std::cref(v), Sequence::dna_alphabet, seq.length() );

  BOOST_CHECK_EQUAL(seq,seq2);
}

BOOST_AUTO_TEST_CASE( odd_dna_string )
{
  std::string seq("ATGCNT");
  Sequence::pack8::vtype v = Sequence::pack8::dna2vtype(std::cref(seq), Sequence::dna_alphabet);
  std::string seq2 = Sequence::pack8::vtype2dna( std::cref(v), Sequence::dna_alphabet, seq.length() );

  BOOST_CHECK_EQUAL(seq,seq2);
}

BOOST_AUTO_TEST_CASE( even_poly_string )
{
  std::string seq("ATGC0");
  Sequence::pack8::vtype v = Sequence::pack8::dna2vtype(std::cref(seq), Sequence::dna_poly_alphabet);
  std::string seq2 = Sequence::pack8::vtype2dna( std::cref(v), Sequence::dna_poly_alphabet, seq.length() );

  BOOST_CHECK_EQUAL(seq,seq2);
}

BOOST_AUTO_TEST_CASE( odd_poly_string )
{
  std::string seq("ATGC01");
  Sequence::pack8::vtype v = Sequence::pack8::dna2vtype(std::cref(seq), Sequence::dna_poly_alphabet);
  std::string seq2 = Sequence::pack8::vtype2dna( std::cref(v), Sequence::dna_poly_alphabet, seq.length() );

  BOOST_CHECK_EQUAL(seq,seq2);
}
