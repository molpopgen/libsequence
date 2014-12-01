#define BOOST_TEST_MODULE FastaConstructors
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <fstream>
#include <boost/test/unit_test.hpp>

std::string name("seqname"),seq("AGCGTAGACAGTAGAGTGAT");

BOOST_AUTO_TEST_CASE( ostream_test )
{
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o("testout.fasta");
  o << f << std::endl;
  o.close();
  std::ifstream in("testout.fasta");
  in >> f2 >> std::ws;
  BOOST_REQUIRE(f==f2);
}

//EOF
