#ifndef __SEQUENCE_IOHELP_HPP__
#define __SEQUENCE_IOHELP_HPP__

#include <string>
#include <zlib.h>

namespace Sequence
{
  namespace IOhelp
  {
    /*!
      Reads from a gzipped stream of ASCI data until whitespace is encountered.
      \return The return value of gzread(gzfile,&ch,sizeof(char))
      \note The buffer is not cleared/emptied by this routine.
    */
    int gzread2ws( gzFile gzfile, std::string & buffer );

    /*!
      Chews through white space is a gzipped stream of ASCI data.
      \return The return value of gzread(gzfile,&ch,sizeof(char))
     */
    int gzreadws( gzFile gzfile );

    int gzreaduntil( gzFile file, const char & until );
  }
}

#endif
