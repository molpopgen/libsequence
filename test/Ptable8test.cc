/*! 
  \file alphabets.cc @brief Unit tests for Sequence::Ptable8
 */
#define BOOST_TEST_MODULE Ptable8test
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Ptable8.hpp>
#include <Sequence/Ptable.hpp>
#include <Sequence/PolySites.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>

//convert, convert back.  All good?
BOOST_AUTO_TEST_CASE( create_1 )
{
  using psite = Sequence::polymorphicSite;
  Sequence::Ptable8 t = { psite(1.,"AAGC"),
			  psite(2.,"ACCA") }; 
  std::string x (  t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}

BOOST_AUTO_TEST_CASE( create_2 )
{
  using psite = Sequence::polymorphicSite;
  std::vector<psite> ps = { psite(1.,"AAGC"),
			   psite(2.,"ACCA") }; 
  Sequence::Ptable8 t(ps);
 
  BOOST_CHECK( ps.empty() );

  std::string x ( t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}

BOOST_AUTO_TEST_CASE( create_3 )
{
  using psite = Sequence::polymorphicSite;
  Sequence::Ptable __t = { psite(1.,"AAGC"),
			   psite(2.,"ACCA") };

  Sequence::Ptable8 t(__t);

  BOOST_CHECK( __t.empty() );

  std::string x ( t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}


BOOST_AUTO_TEST_CASE( create_4 )
{
  using psite = Sequence::polymorphicSite;
  Sequence::Ptable __t = { psite(1.,"AAGC"),
			   psite(2.,"ACCA") };

  Sequence::Ptable8 t(std::cref(__t));

 BOOST_CHECK( !__t.empty() );

 std::string x ( t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}
