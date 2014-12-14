//\file FastaOperations.cc
#define BOOST_TEST_MODULE FastaOperations
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>

//A generic revcom routine written for this test
std::string rcom( const std::string & s )
{
  std::string rv(s);
  std::reverse(rv.begin(),rv.end());
  std::transform(rv.begin(),rv.end(),
		 rv.begin(),
		 [](const char & ch)
		 {
		   switch(ch)
		     {
		     case 'A':
		       return 'T';
		       break;
		     case 'G':
		       return 'C';
		       break;
		     case 'C':
		       return 'G';
		       break;
		     case 'T':
		       return 'A';
		       break;
		     }
		   return 'N';
		 });
  return rv;
}


BOOST_AUTO_TEST_CASE( revcom )
{
  std::string name("seqname"),seq("AGCGTAGACAGTAGAGTGAT");
  Sequence::Fasta f(name,seq);

  Sequence::Fasta f2 = f;
  f2.Revcom();

  BOOST_REQUIRE( f2.second == rcom(seq) );
}

BOOST_AUTO_TEST_CASE( subseq )
{
  std::string name("seqname"),seq("AGCGTAGACAGTAGAGTGAT");
  Sequence::Fasta f(name,seq);

  Sequence::Fasta f3(f);
  f3.Subseq(1,3);

  BOOST_REQUIRE( f3.second == "GCG" );

  f3.Complement();

  BOOST_REQUIRE( f3.second == "CGC" );

  BOOST_REQUIRE( std::string(f3) == "CGC" ); //operator string()

}


BOOST_AUTO_TEST_CASE( gapped )
{
  Sequence::Fasta f3("seqname","GCG");

  BOOST_REQUIRE( !f3.IsGapped() );

  f3.second += '-';

  BOOST_REQUIRE( f3.IsGapped() );

  BOOST_REQUIRE( f3.length() == 4 );

  BOOST_REQUIRE( f3.UngappedLength() == 3 );

  //Remove the gap
  f3.second.erase( f3.second.find('-'), 1 );

  BOOST_REQUIRE( f3.length() == 3 );

  BOOST_REQUIRE( f3.UngappedLength() == 3 );
}

BOOST_AUTO_TEST_CASE( cpp11access_1 )
{
  Sequence::Fasta f3("seqname","GCG");
  for( auto & d  : f3 )
    {
      d = 'A';
    }
  BOOST_REQUIRE_EQUAL(f3.second,"AAA");
}

//EOF
