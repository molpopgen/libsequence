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

    struct gzfillbuffer
    {
      template<typename policy>
      std::pair<std::string,int> operator()( gzFile gzfile, const policy & p )
      {
	char ch;
	std::string rv;
	int gzrv;
	while( (gzrv = gzread(gzfile,&ch,sizeof(char))) != 0 )
	  {
	    if(p(ch))
	      {
		return std::make_pair(rv,gzrv);
	      }
	    else
	      {
		rv += ch;
	      }
	  }
	return std::make_pair(rv,gzrv);
      }
    };

    struct gzreaduntil
    {
      /*!
	Reads data char-by-char until the policy is satisfied or EOF/error
	are encountered, whichever comes first.
	\param gzfile A gzFile
	\param p A function taking a const char & as an argument and returning bool.
	\return The total number of bytes read from the file and the last value returned by gzread;
      */
      template<typename policy>
      std::pair<int,int> operator()( gzFile gzfile, const policy & p )
      {
	char ch;
	int rv = 0;
	int gzrv;
	while( (gzrv = gzread(gzfile,&ch,sizeof(char))) != 0 )
	  {
	    rv += gzrv;
	    if(p(ch))
		{
		  gzungetc(ch,gzfile);
		  return std::make_pair(rv,gzrv);
		}
	  }
	return std::make_pair(rv,gzrv);
      }
    };


    /*!
      Reads from a gzipped stream of ASCI data until whitespace is encountered.
      \return The return value of gzread(gzfile,&ch,sizeof(char))
      \note The buffer is not cleared/emptied by this routine.
    */
    std::pair<int,int> gzread2ws( gzFile gzfile, std::string & buffer );

    /*!
      Chews through white space is a gzipped stream of ASCI data.
      \return The return value of gzread(gzfile,&ch,sizeof(char))
     */
    std::pair<int,int> gzreadws( gzFile gzfile );


    /*!
      Reads until a character is met, or EOF/error are encountered,
      whichever comes first
      \return The return value of gzread(gzfile,&ch,sizeof(char))
    */
    //std::pair<int,int> gzreaduntil( gzFile file, const char & until );

    /*!
      Reads until newline character is encountered or EOF/error, whichever 
      comes first, and returns what was read in a string.

      An empty string occurs when EOF was hit.
    */
    std::pair<std::string,int> gzreadline(  gzFile gzfile );

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
