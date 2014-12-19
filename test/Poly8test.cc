/*! \file alphabets.cc @brief Unit tests for Sequence/Poly8.hpp and namespace Sequence::poly8*/
#define BOOST_TEST_MODULE poly8
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Poly8.hpp>
#include <Sequence/Fasta.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>

//convert, convert back.  All good?
BOOST_AUTO_TEST_CASE( string_1 )
{
  const std::string x("AGCTN");
  Sequence::poly8::vtype x8 = Sequence::poly8::dna2vtype((x));
  std::string x2 = Sequence::poly8::vtype2dna(x8);
  BOOST_CHECK_EQUAL(x,x2);
}

//convert in non-const context
BOOST_AUTO_TEST_CASE( string_2 )
{
  std::string x("AGCTN"),x_temp(x);
  Sequence::poly8::vtype x8 = Sequence::poly8::dna2vtype(x);
  std::string x2 = Sequence::poly8::vtype2dna(x8);
  BOOST_CHECK( x.empty() );  //x should have been cleared
  BOOST_CHECK_EQUAL(x_temp,x2);
}

//From a sequence
BOOST_AUTO_TEST_CASE( Seq_1 )
{
  const Sequence::Fasta x("name","AGCTN");
  Sequence::poly8::vtype x8 = Sequence::poly8::dna2vtype(x);
  std::string x2 = Sequence::poly8::vtype2dna(x8);
  BOOST_CHECK_EQUAL(x.second,x2);
}

//from a non-const sequence
BOOST_AUTO_TEST_CASE( Seq_2 )
{
  Sequence::Fasta x("name","AGCTN");
  auto x_temp = x.second;
  Sequence::poly8::vtype x8 = Sequence::poly8::dna2vtype(x);
  std::string x2 = Sequence::poly8::vtype2dna(x8);
  BOOST_CHECK(x.second.empty());
  BOOST_CHECK_EQUAL( x2, x_temp );
}
