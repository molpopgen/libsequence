#define BOOST_TEST_MODULE FastaConstructors
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <string>
#include <boost/test/unit_test.hpp>

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

//EOF
