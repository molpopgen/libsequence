//! \file polySiteVectorTest.cc @brief Unit tests for Sequence::polySiteVector
#define BOOST_TEST_MODULE polySiteVectorTest

#include <Sequence/polySiteVector.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <boost/test/included/unit_test.hpp>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <functional>

using psite = Sequence::polymorphicSite;
using Ptable = Sequence::polySiteVector;

BOOST_AUTO_TEST_CASE( ptable_remove_1 )
{

  Ptable t = { psite(1.,"AAGC"),
	       psite(2.,"ACZA") }; //site 2 has a non-DNA character
  
  BOOST_CHECK_EQUAL( t.size(), 2 );
  
  t.erase( std::remove_if( t.begin(),
			t.end(),
			[]( const psite & __p ) {
			     return std::find_if(__p.second.begin(),
						 __p.second.end(),
						 Sequence::invalidPolyChar())
			       != __p.second.end();
			   } ),
	   t.end() );
  BOOST_CHECK_EQUAL( t.size(), 1 );
}

BOOST_AUTO_TEST_CASE( ptable_make_from_polytable )
{
 using psite = Sequence::polymorphicSite;
 Ptable t = { psite(1.,"AAGC"),
	      psite(2.,"ACAA") };
 
 Sequence::PolySites ps(t.begin(),t.end());

 BOOST_REQUIRE( std::distance(t.begin(),t.end()) ==
		std::distance(ps.sbegin(),ps.send()) );
 
  auto t_i = t.begin();
  auto ps_i = ps.sbegin();

  while( t_i < t.end() )
    {
      BOOST_CHECK_EQUAL(t_i->first,ps_i->first);
      BOOST_CHECK_EQUAL(t_i->second,ps_i->second);
      ++t_i;
      ++ps_i;
    }

  Ptable t2(Sequence::make_polySiteVector(ps));

  BOOST_REQUIRE( t == t2 );
}
