//\file fastqConstructors.cc
#define BOOST_TEST_MODULE fastqIO

#include <Sequence/fastq.hpp>
#include <fstream>
#include <boost/test/included/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <iostream>
#include <stdexcept>

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
