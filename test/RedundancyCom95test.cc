#define BOOST_TEST_MODULE RedundancyCom95
#define BOOST_TEST_DYN_LINK 

#include <Sequence/RedundancyCom95.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <boost/test/unit_test.hpp>
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
