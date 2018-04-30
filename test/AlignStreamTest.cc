//! \file AlignStreamTest.cc @brief unit tests for Sequence::ClustalW and Sequence::phylipData

#include <Sequence/Fasta.hpp>
#include <Sequence/Clustalw.hpp>
#include <Sequence/phylipData.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <unistd.h>

BOOST_AUTO_TEST_SUITE(AlignStreamTest)
const char * input_data = "data/CG15644-Z.aln";
const char * phylip_input_data = "data/phylip_input.txt";

using clustal = Sequence::ClustalW<Sequence::Fasta>;
using phylip  = Sequence::phylipData<Sequence::Fasta>;


/*
  First, make sure that these use cases compile,
  To check that the static_assert in AlignStream<T>
  is correct
*/
Sequence::ClustalW< std::pair<std::string,std::string > >  clustal_pairs;
Sequence::phylipData< std::pair<std::string,std::string > >  phylip_pairs;

//Test object construction
BOOST_AUTO_TEST_CASE( from_vector1 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };

  BOOST_REQUIRE_NO_THROW( clustal c(vf);
			  BOOST_CHECK_EQUAL(vf.size(),c.size());
			  );
}

BOOST_AUTO_TEST_CASE( from_vector2 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };

  BOOST_REQUIRE_NO_THROW( phylip c(vf);
			  BOOST_CHECK_EQUAL(vf.size(),c.size());
			  );
}

BOOST_AUTO_TEST_CASE( from_vector3 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  BOOST_REQUIRE_NO_THROW (
			  clustal c(std::move(vf));
			  BOOST_CHECK(vf.empty());
			  BOOST_REQUIRE_EQUAL(c.size(),2);
			  );
}

BOOST_AUTO_TEST_CASE( from_vector4 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  BOOST_REQUIRE_NO_THROW (
			  phylip c(std::move(vf));
			  BOOST_CHECK(vf.empty());
			  BOOST_REQUIRE_EQUAL(c.size(),2);
			  );
}

BOOST_AUTO_TEST_CASE( copy_construct1 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  clustal c(std::move(vf));
  BOOST_REQUIRE_NO_THROW(clustal c2(c);
			 BOOST_REQUIRE_EQUAL(c.size(),c2.size());
			 for( decltype(c.size()) i = 0 ; i < c.size() ; ++i )
			   {
			     BOOST_REQUIRE( c[i]==c2[i] );
			   }
			 );
}

BOOST_AUTO_TEST_CASE( copy_construct2 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  phylip c(std::move(vf));
  BOOST_REQUIRE_NO_THROW(phylip c2(c);
			 BOOST_REQUIRE_EQUAL(c.size(),c2.size());
			 for( decltype(c.size()) i = 0 ; i < c.size() ; ++i )
			   {
			     BOOST_REQUIRE( c[i]==c2[i] );
			   }
			 );
}

//Check that forwarding is working.
BOOST_AUTO_TEST_CASE( move_construct1 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  clustal c(std::move(vf));
  BOOST_CHECK_EQUAL(vf.size(),0);
  clustal c2(std::move(c));
  BOOST_REQUIRE_EQUAL(c2.size(),2);
  BOOST_CHECK_EQUAL(c.size(),0);
}

//Check that forwarding is working
BOOST_AUTO_TEST_CASE( move_construct2 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  phylip c(std::move(vf));
  BOOST_CHECK_EQUAL(vf.size(),0);
  phylip c2(std::move(c));
  BOOST_REQUIRE_EQUAL(c2.size(),2);
  BOOST_CHECK_EQUAL(c.size(),0);
}

//Check that forwarding works b/w types
BOOST_AUTO_TEST_CASE( move_construct3 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  clustal c(std::move(vf));
  BOOST_CHECK_EQUAL(vf.size(),0);
  phylip c2(std::move(c));
  BOOST_REQUIRE_EQUAL(c2.size(),2);
  BOOST_CHECK_EQUAL(c.size(),0);
}

//Opposite direction to previous test
BOOST_AUTO_TEST_CASE( move_construct4 )
{
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  phylip c(std::move(vf));
  BOOST_CHECK_EQUAL(vf.size(),0);
  clustal c2(std::move(c));
  BOOST_REQUIRE_EQUAL(c2.size(),2);
  BOOST_CHECK_EQUAL(c.size(),0);
}

BOOST_AUTO_TEST_CASE( should_throw1 )
{
  //Sequences of unequal length = bad mojo
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATCAA") };

  BOOST_REQUIRE_THROW( clustal c(vf), std::runtime_error );
  BOOST_REQUIRE_THROW( phylip  c(vf), std::runtime_error );
  BOOST_REQUIRE_THROW( clustal c( std::vector< Sequence::Fasta >({ Sequence::Fasta("seq1","ATG"),
	    Sequence::Fasta("seq2","ATCAA") })),
    std::runtime_error );
  BOOST_REQUIRE_THROW( clustal c( std::move(std::vector< Sequence::Fasta >({ Sequence::Fasta("seq1","ATG"),
	      Sequence::Fasta("seq2","ATCAA") }))),
    std::runtime_error );
}

BOOST_AUTO_TEST_CASE( should_throw2 ) 
{
  //The data are good
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATC") };
  clustal c(std::move(vf));

  //The user does domething dumb...
  c[1].second += std::string("ATGC");

  //The data are not valid
  BOOST_REQUIRE_EQUAL( c.IsAlignment(),false );

  //And the user did not check and then tries to copy
  BOOST_REQUIRE_THROW( clustal c2(c), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( should_throw3 )
{
  //Sequences of unequal length = bad mojo
  std::vector< Sequence::Fasta > vf = { Sequence::Fasta("seq1","ATG"),
					Sequence::Fasta("seq2","ATCAA") };
  clustal c;
  BOOST_REQUIRE_THROW( c.assign(vf.begin(),
				vf.end()), std::runtime_error );
  //And make sure it works with c++11-style iterator fxns
  BOOST_REQUIRE_THROW( c.assign(vf.cbegin(),
				vf.cend()), std::runtime_error );
}

//Test IO routines

//Simple reading clustalw
BOOST_AUTO_TEST_CASE( clustalw_in )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  BOOST_REQUIRE_NO_THROW( in >> c >> std::ws; );

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , 10 );
}

//simple reading phylip
BOOST_AUTO_TEST_CASE( phylip_in )
{
  std::ifstream in(phylip_input_data);
  phylip p;
  BOOST_REQUIRE_NO_THROW(in >> p >> std::ws;);
}

BOOST_AUTO_TEST_CASE( clustalw_convert )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  BOOST_REQUIRE_NO_THROW( in >> c >> std::ws; );

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , 10 );

  //Copy construction
  phylip p(c);
  
  BOOST_REQUIRE (!p.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , p.size() );
  for( decltype(c.size()) i = 0 ; i < c.size() ; ++i )
    {
      BOOST_REQUIRE_EQUAL( c[i],p[i] );
    }
}

BOOST_AUTO_TEST_CASE( clustalw_convert2 )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  BOOST_REQUIRE_NO_THROW( in >> c >> std::ws; );

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , 10 );

  //Use assignment operator this time
  phylip p = c;
  
  BOOST_REQUIRE (!p.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , p.size() );
  for( decltype(c.size()) i = 0 ; i < c.size() ; ++i )
    {
      BOOST_REQUIRE( c[i]==p[i] );
    }
}

BOOST_AUTO_TEST_CASE( clustalw_reread )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  in >> c >> std::ws;
  in.close();

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , 10 );

  const char * outfile = "clustalw_out_test.aln";
  std::ofstream out(outfile);
  out << c << '\n';
  out.close();

  in.open(outfile);
  clustal c2;
  in >> c2 >> std::ws;

  BOOST_REQUIRE( !c2.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , c2.size() );

  for( decltype(c.size()) i = 0 ; i < c.size() ; ++i )
    {
      BOOST_REQUIRE(c[i]==c2[i]);
    }
  unlink(outfile);
}

BOOST_AUTO_TEST_CASE( convert_write_read )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  in >> c >> std::ws;
  in.close();
  phylip p = c;
  
  const char * outfile = "phylip_test.txt";
  std::ofstream out(outfile);
  BOOST_REQUIRE_NO_THROW( out << p << '\n'; );
  out.close();

  in.open(outfile);

  phylip p2;
  BOOST_REQUIRE_NO_THROW( in >> p2 >> std::ws );
  BOOST_REQUIRE_EQUAL(p.size() , p2.size());

  for( decltype(p.size()) i = 0 ; i < p.size() ; ++i )
    {
      //Note that Sequence::Seq::operator==
      //only compares the sequences, not the 
      //names, so this must still pass
      BOOST_REQUIRE( p[i] == p2[i] );
    }
  unlink(outfile);
}

BOOST_AUTO_TEST_SUITE_END()
//EOF
