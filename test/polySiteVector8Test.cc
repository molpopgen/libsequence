/*! 
  \file alphabets.cc @brief Unit tests for Sequence::polySiteVector8
 */
#define BOOST_TEST_MODULE polySiteVector8test
#define BOOST_TEST_DYN_LINK 

#include <Sequence/polySiteVector8.hpp>
#include <Sequence/polySiteVector.hpp>
#include <Sequence/PolySites.hpp>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>
#include <iostream>

using Ptable = Sequence::polySiteVector;

//convert, convert back.  All good?
BOOST_AUTO_TEST_CASE( create_1 )
{
  using psite = Sequence::polymorphicSite;
  Sequence::polySiteVector8 t = { psite(1.,"AAGC"),
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
  Sequence::polySiteVector8 t(ps);
 
    std::string x ( t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}

BOOST_AUTO_TEST_CASE( create_3 )
{
  using psite = Sequence::polymorphicSite;
  Ptable __t = { psite(1.,"AAGC"),
		 psite(2.,"ACCA") };

  Sequence::polySiteVector8 temp(__t);
  Sequence::polySiteVector8 t(std::move(temp));

  BOOST_CHECK(temp.empty());
  
  std::string x ( t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}


BOOST_AUTO_TEST_CASE( create_4 )
{
  using psite = Sequence::polymorphicSite;
  Ptable __t = { psite(1.,"AAGC"),
			   psite(2.,"ACCA") };

  Sequence::polySiteVector8 t(std::cref(__t));

 BOOST_CHECK( !__t.empty() );

 std::string x ( t.begin()->second.unpack() );

  BOOST_CHECK_EQUAL( x,"AAGC" );

  x = std::string ( (t.begin()+1)->second.unpack() );

  BOOST_CHECK_EQUAL( x,"ACCA" );
}

//Check an issue brought up un UhapsTest.cc
BOOST_AUTO_TEST_CASE( check_nsam_not_equal_nsites )
{
  std::vector<double> pos = {1,2,3,4,5,6,7};
  std::vector<std::string> data = {"AAAAAAA",
				   "AAGAAAA",
				   "CTGAAGA",
				   "NAACTGA",
				   "NAACTGA",
				   "AAACCCA"};
  Sequence::PolySites ps(pos,data);
  Ptable vps(ps.sbegin(),ps.send());
  Sequence::polySiteVector8 t(vps);
  for(unsigned i = 0 ; i < vps.size() ; ++i )
    {
      BOOST_CHECK_EQUAL(t[i].second.unpack(),vps[i].second);
      BOOST_CHECK_EQUAL(t[i].second.unpack(),(ps.sbegin()+i)->second);
    }
}
