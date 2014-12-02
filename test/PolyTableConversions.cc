#define BOOST_TEST_MODULE PolyTableConversions
#define BOOST_TEST_DYN_LINK 

#include <Sequence/PolySites.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <cstdlib>
#include <iostream>

BOOST_AUTO_TEST_CASE( make_a_table_from_strings )
{
  std::vector< std::string >
    data = { ("AATAG"),
	     ("ATTAC") };

  BOOST_REQUIRE( !data.empty() );
  BOOST_REQUIRE_EQUAL( data.size(), 2 );
  BOOST_REQUIRE_EQUAL( data[0].size(),5 );
  BOOST_REQUIRE_EQUAL( data[0].length(),5 );
  BOOST_REQUIRE( Sequence::Alignment::IsAlignment(data) );
  Sequence::PolySites ps(data);
  BOOST_REQUIRE_EQUAL( ps.numsites(), 2 );
}

BOOST_AUTO_TEST_CASE( make_a_table_from_Fasta )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC") };

  BOOST_REQUIRE( !data.empty() );
  BOOST_REQUIRE_EQUAL( data.size(), 2 );
  BOOST_REQUIRE_EQUAL( data[0].length(),5 );
  BOOST_REQUIRE( Sequence::Alignment::IsAlignment(data) );
  Sequence::PolySites ps(data);
  BOOST_REQUIRE_EQUAL( ps.numsites(), 2 );
}

BOOST_AUTO_TEST_CASE( conversion1 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC") };

  Sequence::PolySites ps(data);

  Sequence::SimpleSNP ps2;
  ps2.assign(ps.sbegin(),ps.send());

  BOOST_REQUIRE( ps == ps2 );

  ps.assign(ps2.sbegin(),ps2.send());

  BOOST_REQUIRE( ps == ps2 );
}

BOOST_AUTO_TEST_CASE( conversion2 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };
  
  Sequence::PolySites ps(data);
  ps.Binary();
  Sequence::SimData sd(ps.sbegin(),ps.send());
  
  BOOST_REQUIRE(ps == sd);
}

BOOST_AUTO_TEST_CASE( conversion3 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };
  
  Sequence::PolySites ps(data);
  Sequence::SimpleSNP ss;
  ss.assign(ps.sbegin(),ps.send());
  BOOST_REQUIRE( ss.label(0) != "anc" );
  ss.set_outgroup(true);

  BOOST_REQUIRE_EQUAL( ss.label(0) , "anc" );

  BOOST_REQUIRE_THROW( ss.label(ss.size()), std::out_of_range );
}

BOOST_AUTO_TEST_CASE( conversion2_throw )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };

  //This must not throw, even though our data are not 0/1.
  //Checking what is read in is left up to the programmer
  Sequence::PolySites ps(data);
  BOOST_REQUIRE_NO_THROW(Sequence::SimData sd(ps.sbegin(),ps.send()));
}
