#define BOOST_TEST_MODULE PolySitesIO

#include <Sequence/PolySites.hpp>
#include <boost/test/included/unit_test.hpp>
#include <sstream>
#include <fstream>
#include <unistd.h>
BOOST_AUTO_TEST_CASE( polysites_io )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};
  
  Sequence::PolySites ps(std::move(pos),std::move(data)),ps2;

  std::ostringstream o;
  o << ps << '\n';
  std::istringstream in(o.str());

  BOOST_REQUIRE_NO_THROW( in >> ps2 >> std::ws );

  BOOST_REQUIRE( ps == ps2 );

  const char * fn = "psitesio.txt";

  std::ofstream of(fn);
  of << ps << '\n';
  of.close();
  std::ifstream inf(fn);
  BOOST_REQUIRE_NO_THROW(inf >> ps2 >> std::ws);
  BOOST_REQUIRE( ps == ps2 );
  unlink(fn);
}
