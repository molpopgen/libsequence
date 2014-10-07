//! Various functions for simplifying IO operations
#include <string>
#include <cctype>
#include <zlib.h>

using namespace std;
namespace Sequence
{
  namespace IOhelp
  {
    int gzread2ws( gzFile gzfile, string & buffer )
    {
      char ch;
      int gzrv;
      while( (gzrv = gzread(gzfile,&ch,sizeof(char))) != 0 )
	{
	  if( isspace(ch) || gzeof(gzfile) || gzrv == -1 ) return gzrv;
	  else {
	    buffer += ch;
	  }
	}
      return gzrv;
    }

    int gzreadws( gzFile gzfile )
    {
      char ch;
      int gzrv;
      while( (gzrv = gzread(gzfile,&ch,sizeof(char))) != 0 )
	{
	  if( !isspace(ch) ) 
	    {
	      gzungetc(ch,gzfile);
	      break;
	    }
	}
      return gzrv;
    }

    int gzreaduntil( gzFile gzfile, const char & until )
    {
      char ch;
      int gzrv;
      while( (gzrv = gzread(gzfile,&ch,sizeof(char))) != 0 )
	{
	  if(ch == until)
	    {
	      gzungetc(ch,gzfile);
	      return gzrv;
	    }
	}
      return gzrv;
    }
  }
}
