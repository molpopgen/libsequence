#define BOOST_TEST_MODULE FastaConstructors
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <string>
#include <iostream>
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( fasta_constructors )
{
  std::string name("seqname"),seq("AGCGTAGACAGTAGAGTGAT");

  Sequence::Fasta f; //empty constructor
  BOOST_REQUIRE(f.first.empty());
  BOOST_REQUIRE(f.second.empty());

  f = Sequence::Fasta(name,seq);
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );

  f = Sequence::Fasta(name.c_str(),seq.c_str());
  BOOST_CHECK( f.first == name );
  BOOST_CHECK( f.second == seq );

  Sequence::Fasta f2(f);
  BOOST_REQUIRE(f == f2);

}

//EOF
