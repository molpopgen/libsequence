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
  BOOST_REQUIRE( c.size() == 10 );
}

BOOST_AUTO_TEST_CASE( clustalw_convert )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  in >> c >> std::ws;

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE( c.size() == 10 );

  //Copy construction
  phylip p(c);
  
  BOOST_REQUIRE (!p.empty() );
  BOOST_REQUIRE( c.size() == p.size() );
  for( decltype(c.size()) i = 0 ; i < c.size() ; ++i )
    {
      BOOST_REQUIRE( c[i]==p[i] );
    }
}

BOOST_AUTO_TEST_CASE( clustalw_convert2 )
{
  std::ifstream in(input_data);
  std::string check;
  clustal c;
  in >> c >> std::ws;

  BOOST_REQUIRE( !c.empty() );
  BOOST_REQUIRE( c.size() == 10 );

  //Use assignment operator this time
  phylip p = c;
  
  BOOST_REQUIRE (!p.empty() );
  BOOST_REQUIRE( c.size() == p.size() );
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
  BOOST_REQUIRE( c.size() == 10 );

  const char * outfile = "clustalw_out_test.aln";
  std::ofstream out(outfile);
  out << c << '\n';
  out.close();

  in.open(outfile);
  clustal c2;
  in >> c2 >> std::ws;

  BOOST_REQUIRE( !c2.empty() );
  BOOST_REQUIRE( c.size() == c2.size() );

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
  out << p << '\n';

  unlink(outfile);
}

//EOF
