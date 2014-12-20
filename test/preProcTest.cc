#define BOOST_TEST_MODULE preProcTest
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Poly8.hpp>
#include <Sequence/preProc.hpp>
#include <Sequence/PolySites.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <Sequence/PolySNP.hpp>

BOOST_AUTO_TEST_CASE( create_1 )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  Sequence::preProc ps_pp(ps); //non-move version!

  BOOST_CHECK_EQUAL( ps_pp.ptable.size(),5 );
  BOOST_CHECK_EQUAL( ps_pp.uhaps.size(),4 );

  for( unsigned i = 0 ; i < ps.size() ; ++i )
    {
      std::cerr << ps[i] << ' ' << Sequence::poly8::vtype2dna( ps_pp.uhaps[i] ) << '\n';
      BOOST_CHECK_EQUAL( ps[i], Sequence::poly8::vtype2dna( ps_pp.uhaps[i] ) );
    }
}

BOOST_AUTO_TEST_CASE( create_2 )
{
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data));

  /*
    The "move" version will be 
    2x as fast, use 1/2 the memory,
    etc., by skipping the creation of 
    an expensive temporary
  */
  Sequence::preProc ps_pp(std::move(ps)); 

  BOOST_CHECK_EQUAL( ps_pp.ptable.size(),5 );
  BOOST_CHECK_EQUAL( ps_pp.uhaps.size(),4 );
  //ps should now be empty:
  BOOST_CHECK_EQUAL( ps.empty(), true );
  BOOST_CHECK_EQUAL( ps.numsites(), 0 );
}

BOOST_AUTO_TEST_CASE( create_3 )
{
  using psite = Sequence::polymorphicSite;
  /* Unique haps are:
     AGA
     TCC
     TCN
     ACA
     TCA 
  */
  Sequence::Ptable p = { psite(1.,"ATTATN"),
			 psite(2.,"GCCCCT"),
			 psite(3.,"ACNAAN") };
  Sequence::PolySites ps(p.begin(),p.end());

  Sequence::preProc p_pp(p);

  BOOST_CHECK_EQUAL(p_pp.uhaps.size(),6);

  for(unsigned i = 0 ; i < ps.size() ; ++i )
    {
      BOOST_CHECK_EQUAL(ps[i],Sequence::poly8::vtype2dna(p_pp.uhaps[i]));
    }
}


BOOST_AUTO_TEST_CASE( create_3b )
{
  Sequence::Ptable pt( std::move(  std::vector<Sequence::polymorphicSite>( 1000000,std::make_pair(1.,
												  std::string(1000,'A')) ) ) );
  unsigned j = 0;
  for( auto & s : pt )
    {
      if ( j < 1000 )
	s.second[j++]='T';
      else break;
    }
  Sequence::preProc pp_pt(pt);
  std::cerr << pp_pt.uhaps.size() << '\n';
}

// BOOST_AUTO_TEST_CASE( create_4 )
// {
//   Sequence::PolySites pt ( std::move(std::vector<double>(1000000,1.)),
// 			   std::move(std::vector<std::string>(1000,std::string(1000000,'A')) ));
//   unsigned j = 0;
//   for( auto & s : pt )
//     {
//       s[j++]='T';
//     }
//   std::cerr << "cerr creating PolySNP: ";
//   Sequence::PolySNP apt(&pt);
//   std::cerr << "done\n";
//   std::cerr << apt.ThetaPi() << '\n';
//   //Sequence::preProc pp_pt(std::move(pt));
//   //  std::cerr << pp_pt.ptable.size() << '\n';
// }
