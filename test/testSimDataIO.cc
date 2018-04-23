#define BOOST_TEST_MODULE testSimDataIO

#include <Sequence/SimDataIO.hpp>
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <iterator>
#include <unistd.h>

const char * infile = "data/single_ms.txt";

BOOST_AUTO_TEST_CASE( SimData_read_stream )
{
  std::ifstream in(infile);
  Sequence::SimData d;

  in >> d >> std::ws; 
  BOOST_REQUIRE_EQUAL(d.size(),100);
  BOOST_REQUIRE_EQUAL(d.numsites(),196);
}

BOOST_AUTO_TEST_CASE( SimData_read_stream2 )
{
  std::ifstream in(infile);
  Sequence::SimData d;
  in >> d >> std::ws; 
  in.close();
  BOOST_REQUIRE_EQUAL(d.size(),100);
  BOOST_REQUIRE_EQUAL(d.numsites(),196);

  const char * outfile = "SimDataTest.txt";
  std::ofstream out(outfile);
  //no new line at end, for extra challenge
  out << d << '\n' << d;
  out.close();
  unsigned count = 0;
  in.open(outfile);
  if(!in) std::exit(1);
  BOOST_REQUIRE_NO_THROW (
			  while(!in.eof())
			  {
			    in >> d >> std::ws;
			    ++count;
			  }
			  );
  BOOST_REQUIRE_EQUAL(count,2);
  in.close();
  unlink(outfile);
}

BOOST_AUTO_TEST_CASE( SimData_read_stream3 )
{

  std::ifstream in(infile);
  Sequence::SimData d;
  in >> d >> std::ws; 
  in.close();
  BOOST_REQUIRE_EQUAL(d.size(),100);
  BOOST_REQUIRE_EQUAL(d.numsites(),196);

  const char * outfile = "SimDataTest.txt";
  std::ofstream out(outfile);
  //no new line at end, for extra challenge
  out << d << '\n' << d;
  out.close();
  unsigned count = 0;
  in.open(outfile);
  if(!in) std::exit(1);
  BOOST_REQUIRE_NO_THROW (
			  std::istream_iterator<Sequence::SimData> i(in);
			  for( ; i != std::istream_iterator<Sequence::SimData>() ; ++i )
			    {
			      BOOST_CHECK( *i == d );
			      ++count;  
			    }
			  );
  BOOST_CHECK_EQUAL(count,2);
  in.close();
  unlink(outfile);
}

BOOST_AUTO_TEST_CASE( SimData_read_FILE )
{
  FILE * fp = fopen(infile,"r");
  Sequence::SimData d;
  d.fromfile(fp);
  fclose(fp);

  BOOST_REQUIRE_EQUAL(d.size(),100);
  BOOST_REQUIRE_EQUAL(d.numsites(),196);
}

BOOST_AUTO_TEST_CASE( SimData_read_same )
{
  std::ifstream in(infile);
  Sequence::SimData d,d2;
  in >> d >> std::ws;
  in.close();
  BOOST_REQUIRE_EQUAL(d.size(),100);
  BOOST_REQUIRE_EQUAL(d.numsites(),196);

  FILE * fp = fopen(infile,"r");
  d2.fromfile(fp);
  fclose(fp);

  BOOST_REQUIRE(d == d2);
}

BOOST_AUTO_TEST_CASE( SimData_binaryIO )
{
  std::ifstream in(infile);
  Sequence::SimData d,d2;
  in >> d >> std::ws;
  in.close();

  std::ostringstream out;
  Sequence::write_SimData_binary(out,d);
  std::istringstream in2(out.str());
  d2 = Sequence::read_SimData_binary(in2);

  BOOST_REQUIRE( d==d2 );
}

BOOST_AUTO_TEST_CASE( SimData_gzipIO )
{
  std::ifstream in(infile);
  Sequence::SimData d,d2;
  in >> d >> std::ws;
  in.close();

  const char * outfile = "SimDataTest.gz";
  //gzipped ASCI
  gzFile gzfile = gzopen(outfile,"w");
  auto x = write_SimData_gz(gzfile,d,false);
  gzclose(gzfile);
  BOOST_REQUIRE( x > 0 );
  
  gzfile = gzopen(outfile,"r");
  d2 = Sequence::read_SimData_gz(gzfile);
  gzclose(gzfile);
  BOOST_REQUIRE(d2 == d);

  //gzipped binary
  gzfile = gzopen(outfile,"w");
  x = write_SimData_gz(gzfile,d,true);
  gzclose(gzfile);
  BOOST_REQUIRE( x > 0 );

  gzfile = gzopen(outfile,"r");
  d2 = Sequence::read_SimData_gz(gzfile,true);
  gzclose(gzfile);
  BOOST_REQUIRE(d2 == d);

  unlink(outfile);
}
