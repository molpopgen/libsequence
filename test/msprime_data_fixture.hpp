#ifndef LIBSEQUENCE_TEST_MSPRIME_DATA_FIXTURE_HPP
#define LIBSEQUENCE_TEST_MSPRIME_DATA_FIXTURE_HPP

#include <sstream>
#include "msformatdata.hpp"
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/variant_matrix/msformat.hpp>

struct vmatrix_from_msprime
{
    Sequence::VariantMatrix m;
    Sequence::AlleleCountMatrix c;

    static Sequence::VariantMatrix
    read()
    {
        std::istringstream in(get_msformat_data());
        return Sequence::from_msformat(in);
    }

    vmatrix_from_msprime() : m(read()), c(m) {}
};

struct msprime_stream
{
    std::istringstream in;
    msprime_stream() : in(get_msformat_stream()) {}
};

#endif
