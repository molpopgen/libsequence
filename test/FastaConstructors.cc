//!\ file FastaConstructors.cc

#include <Sequence/Fasta.hpp>
#include <string>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <functional>

struct fasta_constructors_fixture
{
    std::string name, seq;
    fasta_constructors_fixture()
        : name{ "seqname" }, seq{ "AGCGTAGACAGTAGAGTGAT" }
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(FastaConstructorsTest, fasta_constructors_fixture)

BOOST_AUTO_TEST_CASE(empty)
{
    Sequence::Fasta f;
    BOOST_REQUIRE(f.name.empty());
    BOOST_REQUIRE(f.seq.empty());
}

BOOST_AUTO_TEST_CASE(string_con)
{
    Sequence::Fasta f = Sequence::Fasta(name, seq);
    BOOST_CHECK(f.name == name);
    BOOST_CHECK(f.seq == seq);
}

BOOST_AUTO_TEST_CASE(copy_con)
{
    Sequence::Fasta f = Sequence::Fasta(name.c_str(), seq.c_str());
    BOOST_CHECK(f.name == name);
    BOOST_CHECK(f.seq == seq);

    Sequence::Fasta f2(f);
    BOOST_REQUIRE(f == f2);
}

BOOST_AUTO_TEST_CASE(move_con)
{
    Sequence::Fasta f = Sequence::Fasta(name.c_str(), seq.c_str());
    BOOST_CHECK(f.name == name);
    BOOST_CHECK(f.seq == seq);

    Sequence::Fasta f2(std::move(f));
    BOOST_CHECK(f2.name == name);
    BOOST_CHECK(f2.seq == seq);
    BOOST_CHECK(f.length() == 0);
    BOOST_CHECK(f.name.empty());
}

BOOST_AUTO_TEST_CASE(move_con2)
//This "should" work???
{
    std::string a(name), b(seq);
    Sequence::Fasta f = Sequence::Fasta(std::move(a), std::move(b));
    BOOST_CHECK(f.name == name);
    BOOST_CHECK(f.seq == seq);
    BOOST_CHECK(a.empty());
    BOOST_CHECK(b.empty());
}

BOOST_AUTO_TEST_CASE(move_assign)
{
    Sequence::Fasta f = Sequence::Fasta(name, seq);
    BOOST_CHECK(f.name == name);
    BOOST_CHECK(f.seq == seq);

    Sequence::Fasta f2;
    f2 = std::move(f);
    BOOST_CHECK(f2.name == name);
    BOOST_CHECK(f2.seq == seq);
    BOOST_CHECK(f.length() == 0);
    BOOST_CHECK(f.name.empty());
}
BOOST_AUTO_TEST_SUITE_END()
//EOF
