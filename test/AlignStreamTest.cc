#define BOOST_TEST_MODULE AlignStream
#define BOOST_TEST_DYN_LINK

#include <Sequence/Fasta.hpp>
#include <Sequence/Clustalw.hpp>
#include <Sequence/phylipData.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <unistd.h>

const char * input_data = "data/CG15644-Z.aln";

using clustal = Sequence::ClustalW<Sequence::Fasta>;
using phylip  = Sequence::phylipData<Sequence::Fasta>;

//Simple reading
BOOST_AUTO_TEST_CASE( clustalw_in )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  in >> c >> std::ws;

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE_EQUAL( c.size() , 10 );
}

BOOST_AUTO_TEST_CASE( clustalw_convert )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  in >> c >> std::ws;

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
  in >> c >> std::ws;

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
  phylip p = c;
  
  const char * outfile = "phylip_test.txt";
  std::ofstream out(outfile);
  BOOST_REQUIRE_NO_THROW( out << p << '\n'; );
  out.close();

  in.open(outfile);
  phylip p2;
  BOOST_REQUIRE_NO_THROW( in >> p2 >> std::ws );

  BOOST_REQUIRE(!p2.empty());
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

//EOF
