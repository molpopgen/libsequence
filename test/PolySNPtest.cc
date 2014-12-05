#define BOOST_TEST_MODULE PolySNPtest
#define BOOST_TEST_DYN_LINK 

#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <limits>
#include <cmath>

BOOST_AUTO_TEST_CASE( check_empty_table )
{
  Sequence::PolySites d;   //empty data table
  Sequence::PolySNP ad(&d);

  //Estimators of theta = 4Neu, numbers of various types of polymorhpisms
  BOOST_CHECK_EQUAL( ad.ThetaPi(), 0. );
  BOOST_CHECK_EQUAL( ad.ThetaH(), 0. );
  BOOST_CHECK_EQUAL( ad.ThetaL(), 0. );
  BOOST_CHECK_EQUAL( ad.NumPoly(), 0 );
  BOOST_CHECK_EQUAL( ad.NumMutations(), 0 );
  BOOST_CHECK_EQUAL( ad.NumSingletons(), 0 );
  
  //This behavior differs from PolySIM.
  //For PolySIM, if there are no data, the entire block is assumed to be 
  //a matrix of n x L 0s (zeros).  In other words, the entire region is
  //monomorphic for the ancestral state
  //Here, the n x L region has no variation, but you have no idea that
  //the sites are monomorphic for ancestral vs derived character states,
  //hence we return an "impossible" value
  BOOST_CHECK_EQUAL( ad.NumExternalMutations(), std::numeric_limits<unsigned>::max() );

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
