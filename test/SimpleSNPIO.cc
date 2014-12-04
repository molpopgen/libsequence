#define BOOST_TEST_MODULE PolySitesIO
#define BOOST_TEST_DYN_LINK 

#include <Sequence/PolySites.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <boost/test/unit_test.hpp>
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
  
  Sequence::PolySites temp(std::move(pos),std::move(data));

  Sequence::SimpleSNP ps,ps2,ps3;
  ps.assign(temp.sbegin(),temp.send());
  ps3.assign(temp.sbegin(),temp.send());

  std::ostringstream o;
  o << ps << '\n';
  std::istringstream in(o.str());

  BOOST_REQUIRE_NO_THROW( in >> ps2 >> std::ws );

  BOOST_REQUIRE( ps == ps2 );

  const char * fn = "simplesnpio.txt";

  std::ofstream of(fn);
  of << ps << '\n';
  of.close();
  std::ifstream inf(fn);
  BOOST_REQUIRE_NO_THROW(inf >> ps2 >> std::ws);
  BOOST_REQUIRE( ps == ps2 );
  inf.close();
  unlink(fn);

  const char * fn2 = "simplesnpio2.txt";
  //Now, change the outgroup
  ps.set_outgroup(true);
  of.open(fn2);
  of << ps << '\n';
  of.close();
  inf.open(fn2);
  BOOST_REQUIRE_NO_THROW(inf >> ps2 >> std::ws);
  inf.close();
  BOOST_REQUIRE( ps == ps2 );
  BOOST_REQUIRE( ps == ps3 );

  unlink(fn2);
}
