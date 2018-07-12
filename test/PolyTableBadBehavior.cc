/*! \file PolyTableBadBehavior.cc \brief Unit test for really bad things
  
  Most of what is done here are examples of what NOT to do!!!!

  Covers stuff in the following sections of the tutorial:
  - \ref polytable_csi
  - \ref polytable_idiot
 */

#include <Sequence/PolySites.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/polySiteVector.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <functional>

BOOST_AUTO_TEST_SUITE(PolyTableBadBehaviorTest)

BOOST_AUTO_TEST_CASE( exception1 )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAAC"}; //Sequence not same length as pos.size()
  BOOST_REQUIRE_THROW(Sequence::PolySites ps(std::move(pos),std::move(data)), std::runtime_error );
}

// This test became moot in 1.8.9
// BOOST_AUTO_TEST_CASE( badness1 )
// {
//   std::vector<double> pos = {1,2,3,4,5};
//   std::vector<std::string> data = {"AAAAA",
// 				   "AAGAA",
// 				   "CTGAA",
// 				   "NAACT"};

//   Sequence::PolySites ps(std::move(pos),std::move(data)),
//     ps2(Sequence::copyPolyTable(ps));

//   //force the calculation of const_site_iterator's underlying data
//   std::for_each( ps.sbegin(),ps.send(),[]( const Sequence::PolySites::const_site_iterator::value_type & __p ) {} );

//   //Now, remove that haplo with missing data
//   //by skipping PolyTable's iterator functions (this is BAD)
//   ps.second.erase( std::remove_if(ps.second.begin(),
// 				  ps.second.end(),
// 				  [](const std::string & __s) {
// 				    return __s.find('N') != std::string::npos;
// 				  }),
// 		   ps.second.end() );

//   /*
//     Crap!
//     The underlying data did not get reset...
//   */
//   BOOST_CHECK_EQUAL( ps.sbegin()->second,"AACN" );

//   //This will reset the data, as it calls PolyTable::begin in a non-const context
//   auto i = ps.begin();
//   BOOST_CHECK_EQUAL( ps.sbegin()->second,"AAC" );
// }

//Moot in 1.8.9
//Do the above, but correctly this time
// BOOST_AUTO_TEST_CASE( badness1_fix1 )
// {
//   std::vector<double> pos = {1,2,3,4,5};
//   std::vector<std::string> data = {"AAAAA",
// 				   "AAGAA",
// 				   "CTGAA",
// 				   "NAACT"};

//   Sequence::PolySites ps(std::move(pos),std::move(data)),
//     ps2(Sequence::copyPolyTable(ps));

//   //force the calculation of const_site_iterator's underlying data
//   std::for_each( ps.sbegin(),ps.send(),[]( const Sequence::PolySites::const_site_iterator::value_type & __p ) {} );

//   //call PolyTable::begin and end
//   ps.second.erase( std::remove_if(ps.begin(),
// 				  ps.end(),
// 				  [](const std::string & __s) {
// 				    return __s.find('N') != std::string::npos;
// 				  }),
// 		   ps.end() );
//   BOOST_CHECK_EQUAL( ps.sbegin()->second,"AAC" );
// }

//Moot in 1.8.9
// BOOST_AUTO_TEST_CASE( badness1_fix2 )
// {
//   std::vector<double> pos = {1,2,3,4,5};
//   std::vector<std::string> data = {"AAAAA",
// 				   "AAGAA",
// 				   "CTGAA",
// 				   "NAACT"};

//   Sequence::PolySites ps(std::move(pos),std::move(data)),
//     ps2(Sequence::copyPolyTable(ps));

//   //force the calculation of const_site_iterator's underlying data
//   std::for_each( ps.sbegin(),ps.send(),[]( const Sequence::PolySites::const_site_iterator::value_type & __p ) {} );

//   //uh-oh, remove_if is not using PolyTable's functions!
//   ps.second.erase( std::remove_if(ps.second.begin(),
// 				  ps.second.end(),
// 				  [](const std::string & __s) {
// 				    return __s.find('N') != std::string::npos;
// 				  }),
// 		   ps.end() ); //but, we call PolyTable::end here, and the context is non-const, so the next call to sbegin will be fine.
//   BOOST_CHECK_EQUAL( ps.sbegin()->second,"AAC" );
// }

BOOST_AUTO_TEST_CASE( badness2 )
{
  //std::vector<double> pos = {1,2,3,4,5};
  /*
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};
  */
  std::vector<double> pos(1000,1.);
  std::vector<std::string> data(1000,std::string(1000,'A'));
  Sequence::PolySites ps(std::move(pos),std::move(data)),ps2;

  //You now have a table with unequal haplotype lengths
  //This is bad, and is user error.
  ps[0] = std::string("A");

  //Now, your internal data fail the check that all elements
  //are the same length, which may have unintended consequences
  //for later actions.
  BOOST_CHECK_EQUAL(Sequence::Alignment::IsAlignment(std::vector<std::string>(ps.begin(),ps.end())), false);
}
BOOST_AUTO_TEST_SUITE_END()
