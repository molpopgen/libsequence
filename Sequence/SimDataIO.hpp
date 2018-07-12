//! \file SimDataIO.hpp \brief Various I/O functions for Sequence::SimData
#ifndef __SIMDATA_IO_HPP__
#define __SIMDATA_IO_HPP__

#include <Sequence/SimData.hpp>
#include <iosfwd>
#include <zlib.h>
#include <cstring>
#include <cstdio>
#include <unistd.h>

namespace Sequence
{
  /*! \brief Wite a SimData object to a gzFile
    Write a SimData object to a gzFile object (from the C-language zlib library)
    Returns the amount of data written, or -1 on error.
   */
  long long write_SimData_gz( gzFile & file, const SimData & d, const bool & binary = false)__attribute__ ((deprecated));

  /*! \brief Read a SimData object from a gzFile
    Read a SimData object from a gzFile object (from the C-language zlib library)

    \note Only guaranteed compatible with compressed data written using
    write_SimData_gz.  In other words, if you say ms [params] | gzip > file.gz,
    then you cannot expect to use this function to read in the data, and you 
    should use a boost filtering_istream instead.
   */
  SimData read_SimData_gz( gzFile & file, const bool & binary = false )__attribute__ ((deprecated));

  /*! \brief Write a SimData object in binary format to an ostream.
    Write a SimData object in binary format to an ostream.

    The format of the binary data is:
    nsam (uint32_t)
    nsites (uint32_t)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (uint32_t) = # mutations on haplotype i, followed by
    nsites_i values (uint32_t) representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
  */
  void write_SimData_binary( std::ostream & o, const SimData & d )__attribute__ ((deprecated));

  /*!
    Read a SimData object in binary format from an istream

    The format of the binary data is:
    nsam (unsigned)
    nsites (unsigned)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (unsigned) = # mutations on haplotype i, followed by
    nsites_i unsigned values representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
  */
  SimData read_SimData_binary( std::istream & i )__attribute__ ((deprecated));

  /*! \brief  Write a SimData object in binary format to a C-style file descriptor
    Write a SimData object in binary format to a C-style file descriptor
  
    The format of the binary data is:
    nsam (uint32_t)
    nsites (uint32_t)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (uint32_t) = # mutations on haplotype i, followed by
    nsites_i values (uint32_t) representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
   */
  long int write_SimData_binary( int fd , const SimData & d )__attribute__ ((deprecated));

  /*! \brief Write a SimData object in binary format to a C-style file pointer
    Write a SimData object in binary format to a C-style file pointer

    The format of the binary data is:
    nsam (uint32_t)
    nsites (uint32_t)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (uint32_t) = # mutations on haplotype i, followed by
    nsites_i values (uint32_t) representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
   */
  long int write_SimData_binary( FILE * fp , const SimData & d )__attribute__ ((deprecated));
}

#endif
