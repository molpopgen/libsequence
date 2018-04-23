#define BOOST_TEST_MODULE stateCounterTest

#include <Sequence/stateCounter.hpp>
#include <string>
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <algorithm>
#include <functional>

using namespace std;
using namespace Sequence;

BOOST_AUTO_TEST_CASE( test1 )
{
  string x("AGCTN-");
  auto y = for_each(begin(x),end(x),
		    stateCounter());
  BOOST_CHECK_EQUAL(y.a,1);
  BOOST_CHECK_EQUAL(y.g,1);
  BOOST_CHECK_EQUAL(y.c,1);
  BOOST_CHECK_EQUAL(y.t,1);
  BOOST_CHECK_EQUAL(y.n,1);
  BOOST_CHECK_EQUAL(y.gap,1);
  BOOST_CHECK_EQUAL(y.ndna,0);
}
