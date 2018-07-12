/*
  Validate calculation of summary stats via PolySIM
  by independently recoding them all.  Joy.
 */
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <limits>
#include <cmath>

const char * testfile = "data/single_ms.txt";

BOOST_AUTO_TEST_SUITE(PolySIMTest)

double pi( const Sequence::SimData & d )
{
  unsigned ndiffs=0,ncomps=0;

  for( unsigned i = 0 ; i < d.size() - 1 ; ++i )
    {
      for( unsigned j = i + 1 ; j < d.size() ; ++j )
	{
	  ++ncomps;
	  for( unsigned k = 0 ; k < d[i].size() ; ++k )
	    {
	      ndiffs += (d[i][k] != d[j][k]) ? 1 : 0;
	    }
	}
    }
  return double(ndiffs)/double(ncomps);
}

unsigned nsingletons( const Sequence::SimData & d )
{
  if(d.empty())return 0;
  unsigned nsing = 0;
  std::for_each(d.sbegin(),d.send(),
		[&nsing,&d](const Sequence::polymorphicSite & __ps) {
		  auto c = std::count(__ps.second.begin(),__ps.second.end(),'1');
		  nsing += ( unsigned(c) == 1 || unsigned(c) == d.size() - 1 ) ? 1 : 0;
		});
  return nsing;
}

unsigned nexternal( const Sequence::SimData & d )
{
  if(d.empty())return 0;
  unsigned next = 0;
  std::for_each(d.sbegin(),d.send(),
		[&next,&d](const Sequence::polymorphicSite & __ps) {
		  //In "Fu-ese", an external mutation is a derived singleton
		  next += ( std::count(__ps.second.begin(),__ps.second.end(),'1')==1 ) ? 1 : 0;
		});
  return next;
}

//What do we do for no data?
BOOST_AUTO_TEST_CASE( check_empty_table )
{
  Sequence::SimData d;
  Sequence::PolySIM ad(&d);

  //Estimators of theta = 4Neu, numbers of various types of polymorhpisms
  BOOST_CHECK_EQUAL( ad.ThetaPi(), 0. );
  BOOST_CHECK_EQUAL( ad.ThetaH(), 0. );
  BOOST_CHECK_EQUAL( ad.ThetaL(), 0. );
  BOOST_CHECK_EQUAL( ad.NumPoly(), 0 );
  BOOST_CHECK_EQUAL( ad.NumMutations(), 0 );
  BOOST_CHECK_EQUAL( ad.NumSingletons(), 0 );
  BOOST_CHECK_EQUAL( ad.NumExternalMutations(), 0 );

  //Variance on theta estimators that no one in their right mind should use 
  //other than for teaching
  BOOST_CHECK( std::isnan(ad.VarThetaW()) );
  BOOST_CHECK( std::isnan(ad.VarPi()) );
  BOOST_CHECK( std::isnan(ad.SamplingVarPi()) );
  BOOST_CHECK( std::isnan(ad.StochasticVarPi()) );

  //Haplotype-related stats
  /*
    This may seem confusing.  But, and empty PolyTable means no segregating sites.
    But, PolyTables come after filtering.  Further, a SimData is 
    assumed to be based on a sample of size n > 0, thus there is a minimum of 1
    haplotype
  */
  BOOST_CHECK_EQUAL( ad.DandVK(), 1 );  // No. unique haps
  BOOST_CHECK_EQUAL( ad.DandVH(), 0. ); // Haplotype diversity
  BOOST_CHECK( std::isnan(ad.WallsB()) );
  BOOST_CHECK( std::isnan(ad.WallsQ()) );
  BOOST_CHECK_EQUAL( ad.WallsBprime(), 0. );
  BOOST_CHECK_EQUAL( ad.Minrec(), std::numeric_limits<unsigned>::max() );

  //Classic summaries of the site-frequency spectrum
  BOOST_CHECK( std::isnan(ad.TajimasD()) );
  BOOST_CHECK( std::isnan(ad.Dnominator()) );
  BOOST_CHECK( std::isnan(ad.FuLiD()) );
  BOOST_CHECK( std::isnan(ad.FuLiF()) );
  BOOST_CHECK( std::isnan(ad.FuLiDStar()) );
  BOOST_CHECK( std::isnan(ad.FuLiFStar()) );
  BOOST_CHECK( std::isnan(ad.Hprime()) );

  //Other stats
  BOOST_CHECK( std::isnan(ad.HudsonsC()) );
  BOOST_CHECK( ad.Disequilibrium().empty() );
}

BOOST_AUTO_TEST_CASE( pi1 )
{
  Sequence::SimData d;
  std::ifstream in(testfile);
  BOOST_REQUIRE_NO_THROW(in >> d >> std::ws);
  in.close();

  Sequence::PolySIM ad(&d);

  BOOST_REQUIRE_CLOSE( ad.ThetaPi(), pi(d), 0.001 );
  BOOST_REQUIRE_EQUAL( ad.NumPoly(), d.numsites() ); //true by definition of the ms format
}

BOOST_AUTO_TEST_CASE( S1 )
{
  Sequence::SimData d;
  std::ifstream in(testfile);
  BOOST_REQUIRE_NO_THROW(in >> d >> std::ws);
  in.close();

  Sequence::PolySIM ad(&d);

  BOOST_REQUIRE_EQUAL( ad.NumPoly(), d.numsites() ); //true by definition of the ms format
}

BOOST_AUTO_TEST_CASE( singletons1 )
{
  Sequence::SimData d;
  std::ifstream in(testfile);
  BOOST_REQUIRE_NO_THROW(in >> d >> std::ws);
  in.close();

  Sequence::PolySIM ad(&d);

  BOOST_REQUIRE_EQUAL( ad.NumSingletons(), nsingletons(d) ); //true by definition of the ms format
}

BOOST_AUTO_TEST_CASE( external1 )
{
  Sequence::SimData d;
  std::ifstream in(testfile);
  BOOST_REQUIRE_NO_THROW(in >> d >> std::ws);
  in.close();

  Sequence::PolySIM ad(&d);

  BOOST_REQUIRE_EQUAL( ad.NumExternalMutations(), nexternal(d) ); //true by definition of the ms format
}

//number of unique haps
BOOST_AUTO_TEST_CASE( nhaps1 )
{
  Sequence::SimData d;
  std::ifstream in(testfile);
  BOOST_REQUIRE_NO_THROW(in >> d >> std::ws);
  in.close();

  Sequence::PolySIM ad(&d);
  std::set<std::string> uhaps(d.begin(),d.end());
  BOOST_REQUIRE_EQUAL( ad.DandVK(), uhaps.size() );
}

//hap diversity
BOOST_AUTO_TEST_CASE( hapdiv )
{
  Sequence::SimData d;
  std::ifstream in(testfile);
  BOOST_REQUIRE_NO_THROW(in >> d >> std::ws);
  in.close();

  Sequence::PolySIM ad(&d);
  std::set<std::string> uhaps(d.begin(),d.end());
  double hdiv = 0.; //Explicity hap. het
  double hhom = 0.; //Explicity hap. homozygosity.
  for( auto h : uhaps )
    {
      auto c = std::count( d.begin(),d.end(), h );
      hdiv += ( double(c)/double(d.size()) )*(double(d.size()-c)/double(d.size()-1));
      hhom += ( double(c)/double(d.size()) )*(double(c-1)/double(d.size()-1));
    }
  double hdiv2 = 1. - hhom;
  double lseqvalue = ad.DandVH();
  BOOST_REQUIRE_CLOSE( lseqvalue, hdiv, 1e-3);
  BOOST_REQUIRE_CLOSE( lseqvalue, hdiv2, 1e-3);
}
BOOST_AUTO_TEST_SUITE_END()
