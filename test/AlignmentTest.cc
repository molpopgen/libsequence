/*
  Tests for fxns in namespace Sequence::Alignment

  This ns is implmented in terms of T,
  intended to be in thet Sequence::Seq hierarchy.

  Specializations for std::string also exist.

  This module tests bost
*/

#define BOOST_TEST_MODULE AlignmentTest
#define BOOST_TEST_DYN_LINK

#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <boost/test/unit_test.hpp>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <unistd.h>

BOOST_AUTO_TEST_CASE( IsAlignmentFasta )
{
  const char * fn = "GetDataText.txt";

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
			    
			    for(auto n = 0 ; n < vf.size() ; ++n )
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
			    
			    for(auto n = 0 ; n < vf.size() ; ++n )
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
  const char * fn = "GetDataText.txt";

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
