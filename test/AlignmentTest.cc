/*! \file AlignmentTest.cc
  @brief Unit tests for fxns in namespace Sequence::Alignment

  This ns is implmented in terms of T,
  intended to be in thet Sequence::Seq hierarchy.

  Specializations for std::string also exist.

  This module tests both
*/

#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <boost/test/unit_test.hpp>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <unistd.h>

BOOST_AUTO_TEST_SUITE(AlignmentTest)

BOOST_AUTO_TEST_CASE( IsAlignmentFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","AGC") ,
					Sequence::Fasta("seq3","GTC") ,
					Sequence::Fasta("seq4","ATT") ,
					Sequence::Fasta("seq5","AAC") };

  BOOST_REQUIRE_EQUAL( Sequence::Alignment::IsAlignment( vf ), true );
  vf[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::IsAlignment( vf ), true );
  //pull the gaps out in a way that doesn't preserve the alignment
  vf[0].second.erase( std::remove( vf[0].second.begin(),
				   vf[0].second.end(), '-' ), vf[0].second.end() );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::IsAlignment( vf ), false );
}

BOOST_AUTO_TEST_CASE( IsAlignmentString )
{
  std::vector< std::string > vf = { std::string("ATG"),
				    std::string("AGC") ,
				    std::string("GTC") ,
				    std::string("ATT") ,
				    std::string("AAC") };

  BOOST_REQUIRE_EQUAL( Sequence::Alignment::IsAlignment( vf ), true );
  vf[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::IsAlignment( vf ), true );
  //pull the gaps out in a way that doesn't preserve the alignment
  vf[0].erase( std::remove( vf[0].begin(),
			    vf[0].end(), '-' ), vf[0].end() );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::IsAlignment( vf ), false );
}
    
BOOST_AUTO_TEST_CASE( GetDataFasta )
{
  const char * fn = "GetDataText.txt";

  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") },
    vf2;
    
    BOOST_REQUIRE_NO_THROW (
			    std::ofstream o(fn);
			    std::copy(vf.begin(),vf.end(),
				      std::ostream_iterator<Sequence::Fasta>(o,"\n"));
			    o.close();
			    
			    Sequence::Alignment::GetData(vf2,fn);
			    BOOST_REQUIRE(vf == vf2);
			    BOOST_REQUIRE( Sequence::Alignment::IsAlignment(vf2) );
			    unlink(fn);
			    );
}

BOOST_AUTO_TEST_CASE( GetDataFastaStream )
{
  const char * fn = "GetDataText.txt";

  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") },
    vf2;
    
    BOOST_REQUIRE_NO_THROW (
			    std::ofstream o(fn);
			    std::copy(vf.begin(),vf.end(),
				      std::ostream_iterator<Sequence::Fasta>(o,"\n"));
			    o.close();
			    
			    std::ifstream in(fn);
			    Sequence::Alignment::GetData(vf2,in);
			    in.close();
			    BOOST_REQUIRE(vf == vf2);
			    BOOST_REQUIRE( Sequence::Alignment::IsAlignment(vf2) );
			    unlink(fn);
			    );
}

BOOST_AUTO_TEST_CASE( GetDataString )
{
  const char * fn = "GetDataText.txt";

  std::vector< std::string > vf = { "ATG",
				    "ATC" },
    vf2;
    
    BOOST_REQUIRE_NO_THROW (
			    std::ofstream o(fn);
			    std::copy(vf.begin(),vf.end(),
				      std::ostream_iterator<std::string>(o,"\n"));
			    o.close();
			    Sequence::Alignment::GetData(vf2,fn);
			    BOOST_REQUIRE(vf == vf2);
			    BOOST_REQUIRE( Sequence::Alignment::IsAlignment(vf2) );
			    unlink(fn);
			    );
}

BOOST_AUTO_TEST_CASE( GetDataStringStream )
{
  const char * fn = "GetDataText.txt";

  std::vector< std::string > vf = { "ATG",
				    "ATC" },
    vf2;
    
    BOOST_REQUIRE_NO_THROW (
			    std::ofstream o(fn);
			    std::copy(vf.begin(),vf.end(),
				      std::ostream_iterator<std::string>(o,"\n"));
			    o.close();
			    std::ifstream in(fn);
			    Sequence::Alignment::GetData(vf2,in);
			    in.close();
			    BOOST_REQUIRE(vf == vf2);
			    BOOST_REQUIRE( Sequence::Alignment::IsAlignment(vf2) );
			    unlink(fn);
			    );
}

BOOST_AUTO_TEST_CASE( ReadNFasta )
{
  const char * fn = "GetDataText.txt";

  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","AGC") ,
					Sequence::Fasta("seq3","GTC") ,
					Sequence::Fasta("seq4","ATT") ,
					Sequence::Fasta("seq5","AAC") },
    vf2;
    
    BOOST_REQUIRE_NO_THROW (
			    std::ofstream o(fn);
			    std::copy(vf.begin(),vf.end(),
				      std::ostream_iterator<Sequence::Fasta>(o,"\n"));
			    o.close();
			    
			    for(decltype(vf.size()) n = 0 ; n < vf.size() ; ++n )
			      {
				vf2.clear();
				std::ifstream in(fn);
				Sequence::Alignment::ReadNObjects(vf2,n,in);
				BOOST_REQUIRE_EQUAL(vf2.size(),n);
			      }
			    unlink(fn);
			    );
    
    //What happens if you try to read too much in?
    std::ofstream o(fn);
    std::copy(vf.begin(),vf.end(),
	      std::ostream_iterator<Sequence::Fasta>(o,"\n"));
    o.close();
    vf2.clear();
    std::ifstream in(fn);
    BOOST_CHECK_NO_THROW(Sequence::Alignment::ReadNObjects(vf2,vf.size()+1,in));
    BOOST_CHECK_EQUAL(vf2.size(),vf.size());
    in.close();
    unlink(fn);
}

BOOST_AUTO_TEST_CASE( ReadNString )
{
  const char * fn = "GetDataText.txt";

  std::vector< std::string > vf = { std::string("ATG"),
				    std::string("AGC") ,
				    std::string("GTC") ,
				    std::string("ATT") ,
				    std::string("AAC") },
    vf2;
    
    BOOST_REQUIRE_NO_THROW (
			    std::ofstream o(fn);
			    std::copy(vf.begin(),vf.end(),
				      std::ostream_iterator<std::string>(o,"\n"));
			    o.close();
			    
			    for(decltype(vf.size()) n = 0 ; n < vf.size() ; ++n )
			      {
				vf2.clear();
				std::ifstream in(fn);
				Sequence::Alignment::ReadNObjects(vf2,n,in);
				BOOST_REQUIRE_EQUAL(vf2.size(),n);
			      }
			    unlink(fn);
			    );
    
    //What happens if you try to read too much in?
    std::ofstream o(fn);
    std::copy(vf.begin(),vf.end(),
	      std::ostream_iterator<std::string>(o,"\n"));
    o.close();
    vf2.clear();
    std::ifstream in(fn);
    BOOST_CHECK_NO_THROW(Sequence::Alignment::ReadNObjects(vf2,vf.size()+1,in));
    BOOST_CHECK_EQUAL(vf2.size(),vf.size());
    in.close();
    unlink(fn);
}

BOOST_AUTO_TEST_CASE( GappedFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","AGC") ,
					Sequence::Fasta("seq3","GTC") ,
					Sequence::Fasta("seq4","ATT") ,
					Sequence::Fasta("seq5","AAC") };

  BOOST_REQUIRE_EQUAL( Sequence::Alignment::Gapped( vf ), false );
  vf[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::Gapped( vf ), true );
}
 
BOOST_AUTO_TEST_CASE( GappedString )
{
  std::vector< std::string > vf = { std::string("ATG"),
				    std::string("AGC") ,
				    std::string("GTC") ,
				    std::string("ATT") ,
				    std::string("AAC") };
    
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::Gapped( vf ), false );
  vf[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::Gapped( vf ), true );
}

BOOST_AUTO_TEST_CASE( UnGappedLengthFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","AGC") ,
					Sequence::Fasta("seq3","GTC") ,
					Sequence::Fasta("seq4","ATT") ,
					Sequence::Fasta("seq5","AAC") };

  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 3 );
  vf[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 2 );
}

BOOST_AUTO_TEST_CASE( UnGappedLengthString )
{
  std::vector< std::string > vf = { std::string("ATG"),
				    std::string("AGC") ,
				    std::string("GTC") ,
				    std::string("ATT") ,
				    std::string("AAC") };
    
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 3 );
  vf[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ),2 );
}


BOOST_AUTO_TEST_CASE( RemoveGapsFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","AGC") ,
					Sequence::Fasta("seq3","GTC") ,
					Sequence::Fasta("seq4","ATT") ,
					Sequence::Fasta("seq5","AAC") };
  auto vf2 = vf;

  vf2[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );
  Sequence::Alignment::RemoveGaps(vf2);
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), vf2[0].length() );

  for( decltype(vf.size()) i = 0 ; i < vf.size(); ++i )
    {
      //Base 0 should be the same in orig data and copy
      BOOST_REQUIRE_EQUAL( vf[i][0],vf2[i][0] );
      //The gap removal will cause a shift, hence 2 should == 1.
      BOOST_REQUIRE_EQUAL( vf[i][2],vf2[i][1] );
    }
}

BOOST_AUTO_TEST_CASE( RemoveGapsString )
{
  std::vector< std::string > vf = { std::string("ATG"),
				    std::string("AGC") ,
				    std::string("GTC") ,
				    std::string("ATT") ,
				    std::string("AAC") };
  auto vf2 = vf;

  vf2[0][1]='-';
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );
  Sequence::Alignment::RemoveGaps(vf2);
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), vf2[0].length() );

   for( decltype(vf.size()) i = 0 ; i < vf.size(); ++i )
     {
       //Base 0 should be the same in orig data and copy
       BOOST_REQUIRE_EQUAL( vf[i][0],vf2[i][0] );
       //The gap removal will cause a shift, hence 2 should == 1.
       BOOST_REQUIRE_EQUAL( vf[i][2],vf2[i][1] );
     }
}

BOOST_AUTO_TEST_CASE( RemoveTerminalGapsFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","A-TG"),
					Sequence::Fasta("seq2","ACGC") ,
					Sequence::Fasta("seq3","GCTC") ,
					Sequence::Fasta("seq4","ACTT") ,
					Sequence::Fasta("seq5","ACA-") };
  auto vf2 = vf;


  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );
  Sequence::Alignment::RemoveTerminalGaps(vf2);
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );

  //now, vf and vf2 must be the same for the first 3 positions:
  for( decltype(vf.size()) i = 0 ; i < vf.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( vf[i].substr(0,3), vf2[i].second );
      //While we're at it: let's test the explicit type case of a Sequence::Seq to std::string
      BOOST_REQUIRE_EQUAL( vf[i].substr(0,3), std::string(vf2[i]) );
    }
}

BOOST_AUTO_TEST_CASE( RemoveTerminalGapString )
{
  std::vector< std::string > vf = { std::string("A-TG"),
				    std::string("ACGC") ,
				    std::string("GCTC") ,
				    std::string("ACTT") ,
				    std::string("ACA-") };
  auto vf2 = vf;


  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );
  Sequence::Alignment::RemoveTerminalGaps(vf2);
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 2 );

  //now, vf and vf2 must be the same for the first 3 positions:
  for( decltype(vf.size()) i = 0 ; i < vf.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( vf[i].substr(0,3), vf2[i] );
    }
}

BOOST_AUTO_TEST_CASE( TrimFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","A-TG"),
					Sequence::Fasta("seq2","ACGC") ,
					Sequence::Fasta("seq3","GCTC") ,
					Sequence::Fasta("seq4","ACTT") ,
					Sequence::Fasta("seq5","ACA-") };
  
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 2 );
  //We will extract sites 2 and 3 from the alignment using Trim:
  auto vf2 = Sequence::Alignment::Trim( vf, std::vector<int>{1,2} );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 1 );
  //now, vf and vf2 must be the same for the middle 2 positions of vf:
  for( decltype(vf.size()) i = 0 ; i < vf.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( vf[i].second.substr(1,2),vf2[i].second );
    }
}

BOOST_AUTO_TEST_CASE( TrimString )
{
  std::vector< std::string > vf = { std::string("A-TG"),
				    std::string("ACGC") ,
				    std::string("GCTC") ,
				    std::string("ACTT") ,
				    std::string("ACA-") };
  
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 2 );
  //We will extract sites 2 and 3 from the alignment using Trim:
  auto vf2 = Sequence::Alignment::Trim( vf, std::vector<int>{1,2} );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 1 );
  //now, vf and vf2 must be the same for the middle 2 positions of vf:
  for( decltype(vf.size()) i = 0 ; i < vf.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( vf[i].substr(1,2),vf2[i] );
    }
}

BOOST_AUTO_TEST_CASE( TrimComplementFasta )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","A-TG"),
					Sequence::Fasta("seq2","ACGC") ,
					Sequence::Fasta("seq3","GCTC") ,
					Sequence::Fasta("seq4","ACTT") ,
					Sequence::Fasta("seq5","ACA-") };
  
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 2 );
  //We will remove sites 2 and 3 from the alignment using Trim:
  auto vf2 = Sequence::Alignment::TrimComplement( vf, std::vector<int>{1,2} );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 1 );
  //now, vf and vf2 must be the same for the first and last positions of vf
  for( decltype(vf.size()) i = 0 ; i < vf.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( vf[i][0],vf2[i][0] );
      BOOST_REQUIRE_EQUAL( vf[i][3],vf2[i][1] );
    }
}

BOOST_AUTO_TEST_CASE( TrimComplementString )
{
  std::vector< std::string > vf = { std::string("A-TG"),
				    std::string("ACGC") ,
				    std::string("GCTC") ,
				    std::string("ACTT") ,
				    std::string("ACA-") };
  
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf ), 2 );
  //We will remove sites 2 and 3 from the alignment using Trim:
  auto vf2 = Sequence::Alignment::TrimComplement( vf, std::vector<int>{1,2} );
  BOOST_REQUIRE_EQUAL( Sequence::Alignment::UnGappedLength( vf2 ), 1 );
  //now, vf and vf2 must be the same for the first and last positions of vf
  for( decltype(vf.size()) i = 0 ; i < vf.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( vf[i][0],vf2[i][0] );
      BOOST_REQUIRE_EQUAL( vf[i][3],vf2[i][1] );
    }
}

BOOST_AUTO_TEST_CASE( TrimComplementSameResults )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","A-TG"),
					Sequence::Fasta("seq2","ACGC") ,
					Sequence::Fasta("seq3","GCTC") ,
					Sequence::Fasta("seq4","ACTT") ,
					Sequence::Fasta("seq5","ACA-") };
  std::vector< std::string > vs;
  //The magic below works b/c if implicit conversion of Sequence::Seq to std::string
  std::copy( vf.begin(),vf.end(),
	     std::back_inserter(vs) );
 
  auto vf2 = Sequence::Alignment::TrimComplement( vf, std::vector<int>{1,2} );
  auto vs2 = Sequence::Alignment::TrimComplement( vs, std::vector<int>{1,2} );

  for( decltype( vf2.size() ) i = 0 ; i < vf2.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( std::string(vf2[i]),vs2[i] );
      BOOST_REQUIRE_EQUAL( vf2[i].second,vs2[i] );
    }
}
BOOST_AUTO_TEST_SUITE_END()
