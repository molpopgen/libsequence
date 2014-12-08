#define BOOST_TEST_MODULE FastaConstructors
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <string>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <functional>

std::string name("seqname"),seq("AGCGTAGACAGTAGAGTGAT");

BOOST_AUTO_TEST_CASE( empty )
{
  Sequence::Fasta f;
  BOOST_REQUIRE(f.first.empty());
  BOOST_REQUIRE(f.second.empty());
}

BOOST_AUTO_TEST_CASE( string_con )
{
  Sequence::Fasta f = Sequence::Fasta(name,seq);
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );
}

BOOST_AUTO_TEST_CASE( copy_con )
{
  Sequence::Fasta f = Sequence::Fasta(name.c_str(),seq.c_str());
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );
  
  Sequence::Fasta f2(f);
  BOOST_REQUIRE(f == f2);
}

BOOST_AUTO_TEST_CASE( move_con )
{
  Sequence::Fasta f = Sequence::Fasta(name.c_str(),seq.c_str());
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );
  
  Sequence::Fasta f2(std::move(f));
  BOOST_CHECK( f2.first == name );
  BOOST_CHECK( f2.second == seq );
  BOOST_CHECK( f.length() == 0 );
  BOOST_CHECK( f.first.empty() );
}

BOOST_AUTO_TEST_CASE( move_con2 )
//This "should" work???
{
  std::string a(name),b(seq);
  Sequence::Fasta f = Sequence::Fasta(std::move(a),std::move(b));
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );
  BOOST_CHECK( a.empty() );
  BOOST_CHECK( b.empty() );
}

BOOST_AUTO_TEST_CASE( move_assign )
{
  Sequence::Fasta f = Sequence::Fasta(name,seq);
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );

  Sequence::Fasta f2;
  f2 = std::move(f);
  BOOST_CHECK( f2.first == name );
  BOOST_CHECK( f2.second == seq );
  BOOST_CHECK( f.length() == 0 );
  BOOST_CHECK( f.first.empty() );
}
//EOF
