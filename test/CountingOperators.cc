#define BOOST_TEST_MODULE CountingOperators

#include <Sequence/CountingOperators.hpp>
#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <functional>
#include <string>

BOOST_AUTO_TEST_CASE( test_counting_operators_map_plus )
{
  using Sequence::operator+;
  std::map<char,unsigned> baseCounts,baseCounts2;
  baseCounts['A'] = 5;
  baseCounts['G'] = 10;
  baseCounts2['A'] = 11;
  baseCounts2['C'] = 17;
  std::map<char,unsigned> baseCounts3 = baseCounts + baseCounts2;

  BOOST_REQUIRE_EQUAL( baseCounts3['A'], 16 );
  BOOST_REQUIRE_EQUAL( baseCounts3['G'], 10 );
  BOOST_REQUIRE_EQUAL( baseCounts3['C'], 17 );
}

BOOST_AUTO_TEST_CASE( test_counting_operators_map_plus_equal )
{
  using Sequence::operator+=;
  std::map<char,unsigned> baseCounts,baseCounts2;
  baseCounts['A'] = 5;
  baseCounts['G'] = 10;
  baseCounts2['A'] = 11;
  baseCounts2['C'] = 17;
  baseCounts += baseCounts2;

  BOOST_REQUIRE_EQUAL( baseCounts['A'], 16 );
  BOOST_REQUIRE_EQUAL( baseCounts['G'], 10 );
  BOOST_REQUIRE_EQUAL( baseCounts['C'], 17 );
}

BOOST_AUTO_TEST_CASE( test_counting_operators_vector_plus )
{
  using Sequence::operator+;
  std::vector< std::pair<char,unsigned> > baseCounts,baseCounts2;
  baseCounts.push_back(std::make_pair('A',5u));
  baseCounts.push_back(std::make_pair('G',10u));
  baseCounts2.push_back(std::make_pair('A',11u));
  baseCounts2.push_back(std::make_pair('C',17u));

  auto baseCounts3 = baseCounts + baseCounts2;

  std::string bases = {'A','G','C'};

  auto i = std::find_if(baseCounts3.cbegin(),
			baseCounts3.cend(),
			[](const std::pair<char,unsigned> & __p) {
			  return __p.first == 'A';
			});
  BOOST_REQUIRE( i != baseCounts3.cend() );
  BOOST_REQUIRE_EQUAL( i->second,16 );

  i = std::find_if(baseCounts3.cbegin(),
		   baseCounts3.cend(),
		   [](const std::pair<char,unsigned> & __p) {
		     return __p.first == 'G';
		   });
  BOOST_REQUIRE( i != baseCounts3.cend() );
  BOOST_REQUIRE_EQUAL( i->second,10 );
  
  i = std::find_if(baseCounts3.cbegin(),
		   baseCounts3.cend(),
		   [](const std::pair<char,unsigned> & __p) {
		     return __p.first == 'C';
		   });
  BOOST_REQUIRE( i != baseCounts3.cend() );
  BOOST_REQUIRE_EQUAL( i->second,17 );

}
