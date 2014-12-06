#define BOOST_TEST_MODULE SeqConverions
#define BOOST_TEST_DYN_LINK 

#include <Sequence/Fasta.hpp>
#include <Sequence/fastq.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( fastq2fasta )
{
  std::ifstream in("data/data.fastq");
  Sequence::fastq fq;
  Sequence::Fasta fa;
  
  in >> fq >> std::ws;

  fa = fq;

  BOOST_CHECK_EQUAL( fq.first , fa.first );
  BOOST_CHECK_EQUAL( fq.second , fa.second );
}

BOOST_AUTO_TEST_CASE( fastq2fasta2 )
{
  std::ifstream in("data/data.fastq");
  Sequence::fastq fq;
  
  in >> fq >> std::ws;

  Sequence::Fasta fa(fq);

  BOOST_CHECK_EQUAL( fq.first , fa.first );
  BOOST_CHECK_EQUAL( fq.second , fa.second );
}

BOOST_AUTO_TEST_CASE( fastq2fasta3 )
{
  std::ifstream in("data/data.fastq");
  Sequence::fastq fq;
  
  in >> fq >> std::ws;

  Sequence::Fasta fa(std::move(fq));

  BOOST_CHECK (fq.length() == 0);
  BOOST_CHECK (fq.first.empty());
}

BOOST_AUTO_TEST_CASE( fasta2fastq_1 )
{
  Sequence::Fasta fa = {"name","ATGC"};
  Sequence::fastq fq = fa;

  BOOST_CHECK_EQUAL( fq.first , fa.first );
  BOOST_CHECK_EQUAL( fq.second , fa.second );
  BOOST_CHECK( fq.quality.empty() );
}

BOOST_AUTO_TEST_CASE( fasta2fastq_2 )
{
  Sequence::Fasta fa = {"name","ATGC"};
  Sequence::fastq fq = std::move(fa);

  BOOST_CHECK( fq.first == "name" );
  BOOST_CHECK( fq.second == "ATGC" );
  BOOST_CHECK( fq.quality.empty() );
}

BOOST_AUTO_TEST_CASE( fasta2fastq_3 )
{
  Sequence::Fasta fa = {"name","ATGC"};
  Sequence::fastq fq(std::move(fa));

  BOOST_CHECK( fq.first == "name" );
  BOOST_CHECK( fq.second == "ATGC" );
  BOOST_CHECK( fq.quality.empty() );
}
