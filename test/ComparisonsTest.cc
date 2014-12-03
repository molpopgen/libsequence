//Tests for Sequence/Comparisons.hpp
#define BOOST_TEST_MODULE ComparisonsTest
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Comparisons.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>

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
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('A','A'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('A','a'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('a','a'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('A','G'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('A','G'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('g','a'), Sequence::Mutations::Ts );

  BOOST_REQUIRE_EQUAL( Sequence::TsTv('C','T'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('C','t'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('c','t'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('C','T'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('C','T'), Sequence::Mutations::Ts );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('T','c'), Sequence::Mutations::Ts );

  BOOST_REQUIRE_EQUAL( Sequence::TsTv('A','T'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('A','C'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('C','A'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('T','A'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('C','a'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('T','a'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('c','A'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('t','A'), Sequence::Mutations::Tv );

  BOOST_REQUIRE_EQUAL( Sequence::TsTv('G','T'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('G','t'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('g','t'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('G','C'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('G','c'), Sequence::Mutations::Tv );
  BOOST_REQUIRE_EQUAL( Sequence::TsTv('T','g'), Sequence::Mutations::Tv );
}
