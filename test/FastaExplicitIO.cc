#define BOOST_TEST_MODULE FastaExplicitIO
#define BOOST_TEST_DYN_LINK 

#include <Sequence/FastaExplicit.hpp>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>

std::string name("seqname is a seq"),seq("AGCGTAGACAGTAGAGTGAT"),seq_left("AGCGTAGAC"),seq_right("AGTAGAGTGAT");

BOOST_AUTO_TEST_CASE( ostream_test )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  o << f << std::endl;
  o.close();
  std::ifstream in(filename);
  in >> f2 >> std::ws;
  BOOST_REQUIRE_EQUAL(f,f2);
  unlink(filename);
}

BOOST_AUTO_TEST_CASE( ostream_test_3recs )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  o << f << '\n'
    << f << '\n'
    << f << '\n';
  o.close();
  std::ifstream in(filename);
  unsigned count = 0;
  while(! in.eof() )
    {
      in >> f2 >> std::ws;
      BOOST_REQUIRE_EQUAL(f,f2);
      ++count;
    }
  BOOST_REQUIRE_EQUAL( count, 3 );
  //in >> f2 >> std::ws;
  
  unlink(filename);
}

//No newline @ end of data
BOOST_AUTO_TEST_CASE( ostream_test_no_newline_end_of_file )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  o << f << '\n'
    << f << '\n'
    << f;
  o.close();
  std::ifstream in(filename);
  unsigned count = 0;
  while(! in.eof() )
    {
      in >> f2 >> std::ws;
      BOOST_REQUIRE_EQUAL(f,f2);
      ++count;
    }
  BOOST_REQUIRE_EQUAL( count, 3 );
  unlink(filename);
}

BOOST_AUTO_TEST_CASE( ostream_test_newline_within_seq )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  for ( unsigned i = 0 ; i < 3 ; ++i )
    {
      o << '>' << name << '\n'
	<< seq_left << '\n'
	<< seq_right << '\n';
    }
  o.close();
  std::ifstream in(filename);
  unsigned count = 0;
  while(! in.eof() )
    {
      in >> f2 >> std::ws;
      BOOST_REQUIRE_EQUAL(f,f2);
      ++count;
    }
  BOOST_REQUIRE_EQUAL( count, 3 );
  unlink(filename);
}

BOOST_AUTO_TEST_CASE( ostream_test_really_bad_input )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  for ( unsigned i = 0 ; i < 3 ; ++i )
    {
      o << '>' << name << '\n'
	<< seq_left << "\n\n\n\n\n" //Bad bad bad
	<< seq_right << '\n';
    }
  o.close();
  std::ifstream in(filename);
  unsigned count = 0;
  while(! in.eof() )
    {
      in >> f2 >> std::ws;
      BOOST_REQUIRE_EQUAL(f,f2);
      ++count;
    }
  BOOST_REQUIRE_EQUAL( count, 3 );
  unlink(filename);
}

BOOST_AUTO_TEST_CASE( ostream_test_really_bad_input_istream_iterator )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
  Sequence::Fasta f(name,seq),f2;
  std::ofstream o(filename);
  for ( unsigned i = 0 ; i < 3 ; ++i )
    {
      o << '>' << name << '\n'
	<< seq_left << "\n\n\n\n\n" //Bad bad bad
	<< seq_right << '\n';
    }
  o.close();
  std::ifstream in(filename);
  unsigned count = 0;
  std::istream_iterator<Sequence::Fasta> i(in);
  for( ; i != std::istream_iterator<Sequence::Fasta>() ; ++i )
    {
      ++count;  
    }
  BOOST_REQUIRE_EQUAL( count, 3 );
  unlink(filename);
}

BOOST_AUTO_TEST_CASE( exception_test )
{
  const char * filename = "fasta_exp_ostream_test_out.fasta";
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
