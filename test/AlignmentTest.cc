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
			    
			    Sequence::Alignment::GetData<Sequence::Fasta>(vf2,fn);
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
			    Sequence::Alignment::GetData<Sequence::Fasta>(vf2,in);
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
			    Sequence::Alignment::GetData<std::string>(vf2,fn);
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
			    Sequence::Alignment::GetData<std::string>(vf2,in);
			    in.close();
			    BOOST_REQUIRE(vf == vf2);
			    BOOST_REQUIRE( Sequence::Alignment::IsAlignment(vf2) );
			    unlink(fn);
			    );
}
