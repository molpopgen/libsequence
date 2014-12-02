#define BOOST_TEST_MODULE FastaIO
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <unistd.h>

std::string name("seqname"),seq("AGCGTAGACAGTAGAGTGAT");

BOOST_AUTO_TEST_CASE( ostream_test )
{
  const char * filename = "fasta_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  o << f << std::endl;
  o.close();
  std::ifstream in(filename);
  in >> f2 >> std::ws;
  BOOST_REQUIRE(f==f2);
  
  unlink(filename);
}

BOOST_AUTO_TEST_CASE( exception_test )
{
  const char * filename = "fasta_ostream_test_out.fasta";
  std::ofstream out(filename);
  //Write a badly-formatted FASTA record (we forgot the >)
  out << name << '\n' << seq;
  out.close();

  std::ifstream in(filename);
  Sequence::Fasta f;
  BOOST_CHECK_THROW( in >> f >> std::ws, Sequence::badFormat );
  unlink(filename);
}

//EOF
