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
  /*!
    Write a SimData object to a gzFile object (from the C-language zlib library)
   */
  long long write_SimData_gz( gzFile & file, const SimData & d);
  /*!
    Read a SimData object from a gzFile object (from the C-language zlib library)

    \note Only guaranteed compatible with compressed data written using
    write_SimData_gz
   */
  SimData read_SimData_gz( gzFile & file );
  /*!
    Write a SimData object in binary format to an ostream.

    The format of the binary data is:
    nsam (unsigned)
    nsites (unsigned)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (unsigned) = # mutations on haplotype i, followed by
    nsites_i unsigned values representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
  */
  void write_SimData_binary( std::ostream & o, const SimData & d );

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
  SimData read_SimData_binary( std::istream & i );

  /*!
    Write a SimData object in binary format to a C-style file descriptor
  
    The format of the binary data is:
    nsam (unsigned)
    nsites (unsigned)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (unsigned) = # mutations on haplotype i, followed by
    nsites_i unsigned values representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
   */
  int write_SimData_binary( int fd , const SimData & d );

  /*!
    Write a SimData object in binary format to a C-style file pointer

    The format of the binary data is:
    nsam (unsigned)
    nsites (unsigned)
    nsites doubles representing the mutation positions
    Then, there are nsam records containing:
    nsites_i (unsigned) = # mutations on haplotype i, followed by
    nsites_i unsigned values representing the indexes (from 0 to nsites-1) where the derived mutations are
    on haplotype i
   */
  int write_SimData_binary( FILE * fp , const SimData & d );
}

#endif
