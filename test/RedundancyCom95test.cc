#define BOOST_TEST_MODULE RedundancyCom95

#include <Sequence/RedundancyCom95.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <iterator>
#include <unistd.h>

//Test that redundancy values sum to 3 for all non-stop codons + universal code
BOOST_AUTO_TEST_CASE( universal_code_1 )
{
  Sequence::RedundancyCom95 rc;

  unsigned ncodons = 0;

  for( unsigned i = 0 ; i < 4 ; ++i )
    {
      for( unsigned j = 0 ; j < 4 ; ++j )
	{
	  for( unsigned k = 0 ; k < 4 ; ++k )
	    {
	      std::string codon( { Sequence::dna_alphabet[i],
		    Sequence::dna_alphabet[j],
		    Sequence::dna_alphabet[k] } );
	      ++ncodons;
	      //Vals must some to 3 for all but the stop codons
	      if( codon != "TGA" && 
		  codon != "TAG" &&
		  codon != "TAA" )
		{
		  BOOST_CHECK_CLOSE( rc.L0_vals(codon) + rc.L2S_vals(codon) +
				     rc.L2V_vals(codon) + rc.L4_vals(codon) , 3.,1e-6 );
		}
	    }
	}
    }
  BOOST_REQUIRE( ncodons == 64 );
}

BOOST_AUTO_TEST_CASE( AAA )
{
  std::string codon("AAA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( AAC )
{
  std::string codon("AAC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( AAT )
{
  std::string codon("AAT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( AAG )
{
  std::string codon("AAG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( ACA )
{
  std::string codon("ACA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( ACC )
{
  std::string codon("ACC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( ACG )
{
  std::string codon("ACG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( ACT )
{
  std::string codon("ACT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( AGA )
{
  std::string codon("AGA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 1. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( AGC )
{
  std::string codon("AGC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon), 1. );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_EQUAL(rc.First2V(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}

BOOST_AUTO_TEST_CASE( AGG )
{
  std::string codon("AGG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_CLOSE(rc.FirstNon(codon), 2./3., 1e-6 );
  BOOST_CHECK_EQUAL(rc.First2S(codon), 0. );
  BOOST_CHECK_CLOSE(rc.First2V(codon), 1./3., 1e-6 );
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon), 0. );
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon), 0. );
  BOOST_CHECK_EQUAL(rc.Third2S(codon), 1. );
  BOOST_CHECK_EQUAL(rc.Third2V(codon), 0. );
}


BOOST_AUTO_TEST_CASE( AGT )
{
  std::string codon("AGT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( ATA )
{
  std::string codon("ATA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),1.);
}

BOOST_AUTO_TEST_CASE( ATC )
{
  std::string codon("ATC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_CLOSE(rc.ThirdNon(codon),1./3.,1e-6);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_CLOSE(rc.Third2S(codon),1./3.,1e-6);
  BOOST_CHECK_CLOSE(rc.Third2V(codon),1./3.,1e-6);
}

BOOST_AUTO_TEST_CASE( ATG )
{
  std::string codon("ATG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( ATT )
{
  std::string codon("ATT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_CLOSE(rc.ThirdNon(codon),1./3.,1e-6);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_CLOSE(rc.Third2S(codon),1./3.,1e-6);
  BOOST_CHECK_CLOSE(rc.Third2V(codon),1./3.,1e-6);
}

BOOST_AUTO_TEST_CASE( CAA )
{
  std::string codon("CAA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CAC )
{
  std::string codon("CAC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CAG )
{
  std::string codon("CAG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CAT )
{
  std::string codon("CAT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CCA )
{
  std::string codon("CCA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CCC )
{
  std::string codon("CCC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CCG )
{
  std::string codon("CCG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CCT )
{
  std::string codon("CCT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CGA )
{
  std::string codon("CGA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.5);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.5);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CGC )
{
  std::string codon("CGC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CGG )
{
  std::string codon("CGG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_CLOSE(rc.FirstNon(codon),2./3.,1e-6);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_CLOSE(rc.First2V(codon),1./3.,1e-6);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CGT )
{
  std::string codon("CGT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CTA )
{
  std::string codon("CTA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CTC )
{
  std::string codon("CTC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CTG )
{
  std::string codon("CTG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( CTT )
{
  std::string codon("CTT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GAA )
{
  std::string codon("GAA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GAC )
{
  std::string codon("GAC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GAG )
{
  std::string codon("GAG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GAT )
{
  std::string codon("GAT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GCA )
{
  std::string codon("GCA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GCC )
{
  std::string codon("GCC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GCG )
{
  std::string codon("GCG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GCT )
{
  std::string codon("GCT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GGA )
{
  std::string codon("GGA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GGC )
{
  std::string codon("GGC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GGG )
{
  std::string codon("GGG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GGT )
{
  std::string codon("GGT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GTA )
{
  std::string codon("GTA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GTC )
{
  std::string codon("GTC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GTG )
{
  std::string codon("GTG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( GTT )
{
  std::string codon("GTT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TAA )
{
  std::string codon("TAA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TAC )
{
  std::string codon("TAC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TAG )
{
  std::string codon("TAG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TAT )
{
  std::string codon("TAT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TCA )
{
  std::string codon("TCA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TCC )
{
  std::string codon("TCC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TCG )
{
  std::string codon("TCG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TCT )
{
  std::string codon("TCT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TGA )
{
  std::string codon("TGA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TGC )
{
  std::string codon("TGC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TGG )
{
  std::string codon("TGG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TGT )
{
  std::string codon("TGT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TTA )
{
  std::string codon("TTA");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TTC )
{
  std::string codon("TTC");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TTG )
{
  std::string codon("TTG");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}

BOOST_AUTO_TEST_CASE( TTT )
{
  std::string codon("TTT");
  Sequence::RedundancyCom95 rc;
  BOOST_CHECK_EQUAL(rc.FirstNon(codon),1.);
  BOOST_CHECK_EQUAL(rc.First2S(codon),0.);
  BOOST_CHECK_EQUAL(rc.First2V(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdNon(codon),0.);
  BOOST_CHECK_EQUAL(rc.ThirdFour(codon),0.);
  BOOST_CHECK_EQUAL(rc.Third2S(codon),1.);
  BOOST_CHECK_EQUAL(rc.Third2V(codon),0.);
}




