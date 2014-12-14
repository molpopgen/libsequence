//\file fastqConstructors.cc
#define BOOST_TEST_MODULE fastqIO
#define BOOST_TEST_DYN_LINK 

#include <Sequence/fastq.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <iostream>


BOOST_AUTO_TEST_CASE( move_construction )
{
  std::ifstream in("data/data.fastq");
  if (!in)
    {
      std::cerr << "Error, couldn't find input file!\n";
      exit(1);
    }
  Sequence::fastq f;

  in >> f >> std::ws;

  Sequence::fastq f2(std::move(f));

  BOOST_CHECK_EQUAL(f.length(),0);
}
