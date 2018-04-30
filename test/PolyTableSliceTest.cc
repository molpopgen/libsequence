//! \file PolyTableSliceTest.cc @brief Tests for Sequence/PolyTableSlice.hpp
#include <boost/test/unit_test.hpp>
#include <Sequence/SimData.hpp>
#include <iostream>
#include <Sequence/PolyTableSlice.hpp>
#include <vector>
#include <utility>
#include <string>

using namespace std;
using namespace Sequence;

BOOST_AUTO_TEST_SUITE(PolyTableSliceTest)

BOOST_AUTO_TEST_CASE( lastwindows1 )
{
  vector<pair<double,string> > data;
  for(double i = 0.05 ; i < 0.9 ; i += 0.01 )
    data.push_back(make_pair(i,string("001000")));

  SimData d(data.begin(),data.end());
  PolyTableSlice<SimData> w(d.sbegin(),d.send(),0.1,0.001,0.,1.);
  unsigned nwindows = unsigned(1./0.001);
  BOOST_REQUIRE_EQUAL(w.size(),nwindows);
}

BOOST_AUTO_TEST_CASE( nwindows1 )
{
  vector<pair<double,string> > data;
  for(double i = 0.05 ; i < 0.9 ; i += 0.01 )
    data.push_back(make_pair(i,string("001000")));
  
  SimData d(data.begin(),data.end());
  PolyTableSlice<SimData> w(d.sbegin(),d.send(),64);
  unsigned ewindows = std::ceil(double(d.numsites())/64);
  BOOST_REQUIRE_EQUAL(w.size(),std::ceil(double(d.numsites())/double(ewindows)));
  for(auto i = w.cbegin();i!=w.cend();++i)
    {
      auto wi = w.get_slice(i);
      BOOST_CHECK( wi.empty() == false );
    }
}

BOOST_AUTO_TEST_CASE( nwindows2 )
{
  //Make 10x as many SNPs
  vector<pair<double,string> > data;
  for(double i = 0.05 ; i < 0.9 ; i += 0.001 )
    data.push_back(make_pair(i,string("001000")));
  
  SimData d(data.begin(),data.end());
  PolyTableSlice<SimData> w(d.sbegin(),d.send(),64);
  unsigned ewindows = std::ceil(double(d.numsites())/64);
  BOOST_REQUIRE_EQUAL(w.size(),std::ceil(double(d.numsites())/double(ewindows)));
  for(auto i = w.cbegin();i!=w.cend();++i)
    {
      auto wi = w.get_slice(i);
      BOOST_CHECK( wi.empty() == false );
    }
}


BOOST_AUTO_TEST_SUITE_END()
