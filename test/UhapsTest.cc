#define BOOST_TEST_MODULE preProcTest
#define BOOST_TEST_DYN_LINK 


#include <Sequence/SummStats/Uhaps.hpp>
#include <Sequence/PolySites.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <Sequence/PolySNP.hpp>

BOOST_AUTO_TEST_CASE( create_from_polysites )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::Uhaps uh(ps);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

BOOST_AUTO_TEST_CASE( create_from_polysites_move )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::Uhaps uh(ps),
    uh2(std::move(ps));

  BOOST_CHECK_EQUAL(ps.empty(),true);
  BOOST_CHECK_EQUAL(ps.numsites(),0);

  /*
  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), uh2[i].unpack() );
    }
  */
}

BOOST_AUTO_TEST_CASE( create_from_polysitevector8)
{
  // std::vector<double> pos = {1,2,3,4,5};
  // std::vector<std::string> data = {"AAAAA",
  // 				   "AAGAA",
  // 				   "CTGAA",
  // 				   "NAACT",
  // 				   "NAACT"};
  std::vector<double> pos = {1,2,3};
  std::vector<std::string> data  = {"AGA",
				    "TCC",
				    "TCN",
				    "ACA",
				    "TCA"};
  Sequence::PolySites ps(std::move(pos),std::move(data));
  Sequence::polySiteVector8 pt8( ps ); 
				
  Sequence::Uhaps uh(pt8);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

