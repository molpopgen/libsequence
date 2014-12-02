#define BOOST_TEST_MODULE PolyTableTweaking
#define BOOST_TEST_DYN_LINK 

#include <Sequence/PolySites.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/PolyTableManip.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <functional>

//Removal of N 
BOOST_AUTO_TEST_CASE( remove_missing_N )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data)),
    ps2(ps),ps3(ps),ps4(ps);

  BOOST_REQUIRE(ps == ps2);

  BOOST_REQUIRE_EQUAL( ps.numsites() , 5 );
  BOOST_REQUIRE_EQUAL( ps.size(), 4 );

  ps.RemoveMissing();

  BOOST_REQUIRE_EQUAL( ps.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps.size(), 4 );

  //Don't remove missing data from the outgroup
  //outgroup is seq w/missing
  ps2.RemoveMissing(true,3); 

  BOOST_REQUIRE_EQUAL( ps2.numsites() , 5 );
  BOOST_REQUIRE_EQUAL( ps2.size(), 4 );

  //Do remove missing data from the outgroup
  //outgroup is seq w/missing data
  ps2.RemoveMissing(false,3);

  BOOST_REQUIRE_EQUAL( ps2.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps2.size(), 4 );

  //Don't remove missing data based on outgroup
  //outgroup does not have missing data
  ps3.RemoveMissing(false,0); 
  BOOST_REQUIRE_EQUAL( ps3.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps3.size(), 4 );

  BOOST_REQUIRE( ps2 == ps3 );

  //Do remove missing data based on outgroup
  //outgroup does not have missing data
  ps4.RemoveMissing(false,0); 
  BOOST_REQUIRE_EQUAL( ps4.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps4.size(), 4 );

  BOOST_REQUIRE( ps3 == ps4 );
}

//What about lower-case data?
BOOST_AUTO_TEST_CASE( remove_missing_n )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"aaaaa",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  for( auto & d : ps )
    {
      std::transform(d.begin(),
		     d.end(),
		     d.begin(),
		     [](char & ch) { return std::tolower(ch); });
    }

  Sequence::PolySites ps2(ps),ps3(ps),ps4(ps);

  BOOST_REQUIRE(ps == ps2);

  BOOST_REQUIRE_EQUAL( ps.numsites() , 5 );
  BOOST_REQUIRE_EQUAL( ps.size(), 4 );

  ps.RemoveMissing();

  BOOST_REQUIRE_EQUAL( ps.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps.size(), 4 );

  //Don't remove missing data from the outgroup
  //outgroup is seq w/missing
  ps2.RemoveMissing(true,3); 

  BOOST_REQUIRE_EQUAL( ps2.numsites() , 5 );
  BOOST_REQUIRE_EQUAL( ps2.size(), 4 );

  //Do remove missing data from the outgroup
  //outgroup is seq w/missing data
  ps2.RemoveMissing(false,3);

  BOOST_REQUIRE_EQUAL( ps2.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps2.size(), 4 );

  //Don't remove missing data based on outgroup
  //outgroup does not have missing data
  ps3.RemoveMissing(false,0); 
  BOOST_REQUIRE_EQUAL( ps3.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps3.size(), 4 );

  BOOST_REQUIRE( ps2 == ps3 );

  //Do remove missing data based on outgroup
  //outgroup does not have missing data
  ps4.RemoveMissing(false,0); 
  BOOST_REQUIRE_EQUAL( ps4.numsites() , 4 );
  BOOST_REQUIRE_EQUAL( ps4.size(), 4 );

  BOOST_REQUIRE( ps3 == ps4 );
}

//Like "nextgen" data, right?
BOOST_AUTO_TEST_CASE( remove_missing_extreme )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"aaaaa",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  //convert it all to missing
  for( auto & d : ps )
    {
      std::transform(d.begin(),
		     d.end(),
		     d.begin(),
		     [](char & ch) { return 'N'; });
    }

  ps.RemoveMissing();
  BOOST_REQUIRE(ps.empty());
  BOOST_REQUIRE_EQUAL(ps.numsites(),0);
  BOOST_REQUIRE_EQUAL(ps.size(),0);
}

//operator== is case-sensitive!
BOOST_AUTO_TEST_CASE( to_lower )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data)),
    ps2(ps);

  for( auto & d : ps )
    {
      std::transform(d.begin(),
		     d.end(),
		     d.begin(),
		     [](char & ch) { return std::tolower(ch); });
    }

  //Now, ps and ps2 will not be equal
  BOOST_REQUIRE( ps != ps2 );

  for( auto & d : ps )
    {
      std::transform(d.begin(),
		     d.end(),
		     d.begin(),
		     [](char & ch) { return std::toupper(ch); });
    }

  BOOST_REQUIRE( ps == ps2 );
}

//Remove sites based on minor allele freq
BOOST_AUTO_TEST_CASE( remove_maf )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  ps.ApplyFreqFilter(2);

  BOOST_REQUIRE_EQUAL( ps.numsites(), 1 );
  BOOST_REQUIRE( !ps.empty() );

  ps.ApplyFreqFilter(3);

  BOOST_REQUIRE_EQUAL( ps.numsites(), 0 );
  BOOST_REQUIRE_EQUAL( ps.size(), 0 );
  BOOST_REQUIRE( ps.empty() );
}
