#include <config.h>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

//From config.h
static const std::string LIBSEQ_VERSION(VERSION);

int main(int argc, char ** argv)
{
  if(argc==1)
    {
      cerr << "usage:\n"
	   << "\t--version\tPrint out version number and exit\n";
      exit(EXIT_SUCCESS);
    }

  string av1(argv[1]);
  if( av1 == "--version" ) cout << LIBSEQ_VERSION << '\n';

  exit(EXIT_SUCCESS);
}
