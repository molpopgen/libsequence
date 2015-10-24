#define BOOST_TEST_MODULE UhapsTest
#define BOOST_TEST_DYN_LINK 


#include <Sequence/SummStats/Uhaps.hpp>
#include <Sequence/PolySites.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <Sequence/PolySNP.hpp>

BOOST_AUTO_TEST_CASE( create_from_polysites_odd )
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

BOOST_AUTO_TEST_CASE( create_from_polysites_even )
{
  std::vector<double> pos = {1,2,3,4,5,6};
  std::vector<std::string> data = {"AAAAAA",
				   "AAGAAA",
				   "CTGAAG",
				   "NAACTG",
				   "NAACTG",
				   "AAACCC"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::Uhaps uh(ps);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

BOOST_AUTO_TEST_CASE( create_from_polysites_nsites_not_equal_nsam )
{
  std::vector<double> pos = {1,2,3};
  std::vector<std::string> data = {"AAA",
				   "AAG",
				   "CTG",
				   "NAA",
				   "NAA",
				   "AAA"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::Uhaps uh(ps);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

BOOST_AUTO_TEST_CASE( create_from_polysitevector8_odd)
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT",
				   "NAACT"};
 
  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::polySiteVector8 pt8( ps ); 
				
  Sequence::Uhaps uh(pt8);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

BOOST_AUTO_TEST_CASE( create_from_polysitevector8_even)
{
  std::vector<double> pos = {1,2,3,4,5,6};
  std::vector<std::string> data = {"AAAAAA",
				   "AAGAAA",
				   "CTGAAG",
				   "NAACTG",
				   "NAACTG",
				   "AAACCC"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::polySiteVector8 pt8( ps ); 
				
  Sequence::Uhaps uh(pt8);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

BOOST_AUTO_TEST_CASE( create_from_polysitevector8_nsites_not_equal_nsam)
{
  std::vector<double> pos = {1,2,3,4,5,6,7};
  std::vector<std::string> data = {"AAAAAAA",
				   "AAGAAAA",
				   "CTGAAGA",
				   "NAACTGA",
				   "NAACTGA",
				   "AAACCCA"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::polySiteVector8 pt8( ps ); 
				
  Sequence::Uhaps uh(pt8);

  for( unsigned i = 0 ; i < uh.size() ; ++i )
    {
      BOOST_CHECK_EQUAL( uh[i].unpack(), ps[i] );
    }
}

