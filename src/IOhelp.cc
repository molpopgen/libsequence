//! Various functions for simplifying IO operations
#include <Sequence/IOhelp.hpp>
#include <cctype>
using namespace std;

namespace Sequence
{
  namespace IOhelp
  {
    pair<int,int> gzread2ws( gzFile gzfile, string & buffer )
    {
      return gzreaduntil()(gzfile,buffer,[](const char & ch){ return std::isspace(ch); } );
    }

    pair<int,int> gzreadws( gzFile gzfile )
    {
      return gzreaduntil()(gzfile,[](const char & ch){ return !std::isspace(ch); } );
    }

    pair<string,int> gzreadline(  gzFile gzfile )
    {
      return gzfillbuffer()(gzfile,[](const char & ch){ return ch == '\n'; });
    }

    template<>
    void writeBin<std::string>(ostream & o, std::string const & s)
    {
      o.write(s.c_str(),s.size());
    }
  }
}
