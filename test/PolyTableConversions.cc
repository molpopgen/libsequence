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

BOOST_AUTO_TEST_CASE( const_construction_1 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };

  Sequence::PolySites ps(data);
  std::vector<double> pos;
  std::vector<std::string> snps;
  std::copy( ps.pbegin(), ps.pend(),
	     std::back_inserter(pos) );
  std::copy( ps.begin(), ps.end(),
	     std::back_inserter(snps) );

  //Is this picking the const constructor,
  //or the move constructor?
  Sequence::PolySites ps2(pos,
			  snps);

  BOOST_CHECK_EQUAL( pos.size(), ps2.numsites() );
  BOOST_CHECK_EQUAL( snps.size(), ps2.size() );
  BOOST_REQUIRE_EQUAL( ps.numsites() , ps2.numsites() );
  BOOST_REQUIRE_EQUAL( ps.size() , ps2.size() );
  BOOST_REQUIRE( ps == ps2 );
}

BOOST_AUTO_TEST_CASE( const_construction_2 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };

  Sequence::PolySites ps(data);
  std::vector<double> pos;
  std::vector<std::string> snps;
  std::copy( ps.pbegin(), ps.pend(),
	     std::back_inserter(pos) );
  std::copy( ps.begin(), ps.end(),
	     std::back_inserter(snps) );

  //These next two represent IQ test failures by a user

  //Is this picking the const constructor,
  //or the move constructor?
  //It should pick the const
  Sequence::PolySites ps2(std::cref(pos),
			  std::move(snps));

  BOOST_CHECK_EQUAL( pos.size(), ps2.numsites() );
  BOOST_CHECK_EQUAL( snps.size(), ps2.size() );
  BOOST_REQUIRE_EQUAL( ps.numsites() , ps2.numsites() );
  BOOST_REQUIRE_EQUAL( ps.size() , ps2.size() );
  BOOST_REQUIRE( ps == ps2 );

  //Test the reverse
  ps2 = Sequence::PolySites(std::move(pos),
			    std::cref(snps));

  BOOST_CHECK_EQUAL( pos.size(), ps2.numsites() );
  BOOST_CHECK_EQUAL( snps.size(), ps2.size() );
  BOOST_REQUIRE_EQUAL( ps.numsites() , ps2.numsites() );
  BOOST_REQUIRE_EQUAL( ps.size() , ps2.size() );
  BOOST_REQUIRE( ps == ps2 );

}

BOOST_AUTO_TEST_CASE( move_construction_1 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };

  Sequence::PolySites ps(data);
  std::vector<double> pos;
  std::vector<std::string> snps;
  std::copy( ps.pbegin(), ps.pend(),
	     std::back_inserter(pos) );
  std::copy( ps.begin(), ps.end(),
	     std::back_inserter(snps) );

  //This should pick the && constructor
  Sequence::PolySites ps2(std::move(pos),
			  std::move(snps));

  BOOST_CHECK_EQUAL( pos.size(), 0 );
  BOOST_CHECK_EQUAL( snps.size(), 0 );
  BOOST_REQUIRE_EQUAL( ps.numsites() , ps2.numsites() );
  BOOST_REQUIRE_EQUAL( ps.size() , ps2.size() );
  BOOST_REQUIRE( ps == ps2 );
}

BOOST_AUTO_TEST_CASE( move_construction_2 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };

  Sequence::PolySites ps(data);

  Sequence::PolySites ps2(std::move(ps));

  BOOST_CHECK(ps.empty());
}

BOOST_AUTO_TEST_CASE( move_construction_3 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };

  Sequence::PolySites ps(data);

  Sequence::PolySites ps2 = std::move(ps);

  BOOST_CHECK(ps.empty());
}

BOOST_AUTO_TEST_CASE( move_construction_4 )
{
  std::ifstream in("data/single_ms.txt");
  Sequence::SimData d;
  in >> d >> std::ws;
  Sequence::SimData d2(std::move(d));
  BOOST_CHECK(d.empty());
}

BOOST_AUTO_TEST_CASE( move_construction_5 )
{
  std::ifstream in("data/single_ms.txt");
  Sequence::SimData d;
  in >> d >> std::ws;
  Sequence::SimData d2 = std::move(d);
  BOOST_CHECK(d.empty());
}

//PolySites can be created via move
BOOST_AUTO_TEST_CASE( move_conversion_1 )
{
  std::ifstream in("data/single_ms.txt");
  Sequence::SimData d;
  in >> d >> std::ws;
  auto nsam = d.size();
  auto nsites = d.numsites();
  Sequence::PolySites d2 = std::move(d);
  BOOST_CHECK(d.empty());
  BOOST_CHECK_EQUAL( d2.size(), nsam );
  BOOST_CHECK_EQUAL( d2.numsites(), nsites );
 }

BOOST_AUTO_TEST_CASE( move_conversion_2 )
{
  std::ifstream in("data/single_ms.txt");
  Sequence::SimData d;
  in >> d >> std::ws;
  auto nsam = d.size();
  auto nsites = d.numsites();
  Sequence::PolySites d2(std::move(d));
  BOOST_CHECK(d.empty());
  BOOST_CHECK_EQUAL( d2.size(), nsam );
  BOOST_CHECK_EQUAL( d2.numsites(), nsites );
 }

//Making SimData via move
BOOST_AUTO_TEST_CASE( move_conversion_3 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };
  
  Sequence::PolySites ps(data);
  ps.Binary();
  Sequence::SimData sd(std::move(ps));

  BOOST_CHECK(ps.empty());
  BOOST_CHECK_EQUAL(sd.size(),3);
  BOOST_CHECK_EQUAL(sd.numsites(),2);
}

BOOST_AUTO_TEST_CASE( move_conversion_4 )
{
  std::vector< Sequence::Fasta >
    data = { Sequence::Fasta("seq1","AATAG"),
	     Sequence::Fasta("seq2","ATTAC"),
	     Sequence::Fasta("seq3","AATAC") };
  
  Sequence::PolySites ps(data);
  ps.Binary();
  Sequence::SimData sd = std::move(ps);

  BOOST_CHECK(ps.empty());
  BOOST_CHECK_EQUAL(sd.size(),3);
  BOOST_CHECK_EQUAL(sd.numsites(),2);
}

BOOST_AUTO_TEST_CASE( assign_conversion_1 )
{
  std::ifstream in("data/single_ms.txt");
  Sequence::SimData d;
  in >> d >> std::ws;
  auto nsam = d.size();
  auto nsites = d.numsites();
  Sequence::PolySites d2 = d;
  BOOST_CHECK(!d.empty());
  BOOST_CHECK(d == d2);
  BOOST_CHECK_EQUAL( d2.size(), nsam );
  BOOST_CHECK_EQUAL( d2.numsites(), nsites );
 }

BOOST_AUTO_TEST_CASE( assign_conversion_2 )
{
  std::ifstream in("data/single_ms.txt");
  Sequence::SimData d;
  in >> d >> std::ws;
  auto nsam = d.size();
  auto nsites = d.numsites();
  Sequence::PolySites d2 = d;
  Sequence::SimData d3 = d2;
  BOOST_CHECK(d == d3);
 }
