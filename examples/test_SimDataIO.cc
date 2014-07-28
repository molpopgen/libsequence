#include <Sequence/SimDataIO.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Sequence;

int main( int argc, char ** argv )
{
  SimData d;

  cin >> d >> ws;

  gzFile gzf = gzopen("test_zlib_out.gz","w");
  write_SimData_gz(gzf, d);
  gzclose(gzf);

  //now, try to read it
  gzf = gzopen("test_zlib_out.gz","r");
  SimData d2 = read_SimData_gz(gzf);
  gzclose(gzf);

  //write it in binary
  ofstream obin("test_binary_out.bin",ios::binary);
  write_SimData_binary(obin,d);
  obin.close();

  //read it
  ifstream ibin("test_binary_out.bin",ios::binary);
  SimData d5 = read_SimData_binary(ibin);
  ibin.close();
  

  cout << d << '\n' << d2 << '\n' << d5 << '\n';
}
