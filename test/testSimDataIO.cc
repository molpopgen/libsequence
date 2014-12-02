#define BOOST_TEST_MODULE FastaIO
#define BOOST_TEST_DYN_LINK 

#include <Sequence/SimDataIO.hpp>
#include <fstream>
#include <cstdio>
#include <boost/test/unit_test.hpp>
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

  BOOST_REQUIRE_EQUAL(d,d2);
}
