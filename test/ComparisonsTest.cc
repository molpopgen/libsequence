//! \file ComparisonsTests.cc @brief Tests for Sequence/Comparisons.hpp
#define BOOST_TEST_MODULE ComparisonsTest

#include <Sequence/Comparisons.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <iterator>

BOOST_AUTO_TEST_CASE( notagap ) //Silly!
{
  BOOST_REQUIRE_EQUAL( Sequence::NotAGap('-'), false );
}

BOOST_AUTO_TEST_CASE( gapped1 )
{
  std::string seq1("ATGATG---ACA");
  
  BOOST_REQUIRE_EQUAL( Sequence::Gapped(seq1), true );

  std::remove(seq1.begin(),seq1.end(),'-');

  BOOST_REQUIRE_EQUAL( Sequence::Gapped(seq1), false );
}

BOOST_AUTO_TEST_CASE( gapped2 )
{
  std::string seq1("ATGATG---ACA");
  
  auto x = seq1.find('-');

  BOOST_CHECK_PREDICATE(std::not_equal_to<decltype(x)>(),(x)(std::string::npos));

  BOOST_REQUIRE_EQUAL( Sequence::Gapped(seq1.begin(),seq1.end()), true );
  BOOST_REQUIRE_EQUAL( Sequence::Gapped(seq1.begin(),seq1.end(),'x'), false );
  BOOST_REQUIRE_EQUAL( Sequence::Gapped(seq1.begin(),seq1.begin()+x), false );

  std::remove(seq1.begin(),seq1.end(),'-');

  BOOST_REQUIRE_EQUAL( Sequence::Gapped(seq1.begin(),seq1.end()), false );
}

BOOST_AUTO_TEST_CASE( numdiffs1 )
{
  std::string seq1("ATGATG---ACA"),
    seq2(seq1);

  BOOST_REQUIRE_EQUAL( Sequence::NumDiffs(seq1,seq2), 0 );

  //Replace gap characters with N
  //Yes, this is done in the most obtuse way possible... :)
  seq2.replace( seq1.find('-'),
		std::count(seq1.begin(),seq1.end(),'-'), 
		std::string(std::count(seq1.begin(),seq1.end(),'-'),'N' ).c_str());

  BOOST_REQUIRE_EQUAL( Sequence::NumDiffs(seq1,seq2), 0 );
  //Allow missing data to be counted as differences
  BOOST_REQUIRE_EQUAL( Sequence::NumDiffs(seq1,seq2,false), 3 );
}

BOOST_AUTO_TEST_CASE( numdiffs2 )
{
  std::string seq1("ATGATG---ACA"),
    seq2(seq1);

  seq2.erase(std::remove(seq2.begin(),seq2.end(),'-'),seq2.end());

  /*
    With the gaps removed, the two sequences are:
    ATGATG---ACA
    ATGATGACA

    NumDiffs must now return -1
  */
  BOOST_REQUIRE_EQUAL( Sequence::NumDiffs(seq1,seq2), -1 );
}

BOOST_AUTO_TEST_CASE( tstv1 )
{
  //The type casting to int is a hack b/c strongly-type enums are not implicitly-convertible to something printable.
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('A','A')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('A','a')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('a','a')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('A','G')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('A','G')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('g','a')), int(Sequence::Mutations::Ts) );

  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('C','T')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('C','t')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('c','t')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('C','T')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('C','T')), int(Sequence::Mutations::Ts) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('T','c')), int(Sequence::Mutations::Ts) );

  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('A','T')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('A','C')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('C','A')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('T','A')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('C','a')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('T','a')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('c','A')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('t','A')), int(Sequence::Mutations::Tv) );

  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('G','T')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('G','t')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('g','t')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('G','C')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('G','c')), int(Sequence::Mutations::Tv) );
  BOOST_REQUIRE_EQUAL( int(Sequence::TsTv('T','g')), int(Sequence::Mutations::Tv) );

  BOOST_REQUIRE_THROW( Sequence::TsTv('G','R'), std::runtime_error );
  BOOST_REQUIRE_THROW( Sequence::TsTv('G','z'), std::runtime_error );
}
