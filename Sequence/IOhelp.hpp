#ifndef __SEQUENCE_IOHELP_HPP__
#define __SEQUENCE_IOHELP_HPP__

#include <string>
#include <istream>
#include <ostream>
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

    /*!
      Reads until a character is met, or EOF/error are encountered,
      whichever comes first
      \return The return value of gzread(gzfile,&ch,sizeof(char))
    */
    int gzreaduntil( gzFile file, const char & until );

    /*!
      Reads from a binary stream by calling
      i.read(reinterpret_cast<char *>(rv),howmany*sizeof(T));

      \note If rv == nullptr, the function exits without reading
     */
    template<typename T>
    void readBin(std::istream & i, T * rv, const size_t & howmany)
    {
      if ( rv == nullptr ) return;
      i.read(reinterpret_cast<char *>(rv),std::streamsize(howmany*sizeof(T)));
    }

    /*!
      Writes to a binary stream by calling
      o.write( reinterpret_cast<const char *>(t), howmany*sizeof(T) );

      \note Returns without readin if t == nullptr
     */
    template<typename T> 
    void writeBin(std::ostream & o, T const * t, const size_t & howmany = 1)
    {
      if(t == nullptr) return;
      o.write( reinterpret_cast<const char *>(t), std::streamsize(howmany*sizeof(T)));
    }
    
    /*!
      Writes to a binary stream by calling
      o.write( reinterpret_cast<const char *>(&t), sizeof(T) );
    */
    template<typename T> 
    void writeBin(std::ostream & o, T const & t)
    {
      o.write( reinterpret_cast<const char *>(&t), sizeof(T) );
    }
    
    //Specialization for std::string
    template<>
    void writeBin<std::string>(std::ostream & o, std::string const & s);
  }
}

#endif
