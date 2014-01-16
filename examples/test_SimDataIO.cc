#include <Sequence/SimDataIO.hpp>
#include <iostream>
#include <fstream>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace Sequence;
using namespace boost::iostreams;

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

  //now, write it via boost
  filtering_ostream o;
  o.push(gzip_compressor());
  o.push(file_sink("test_boostgz_out.gz",ios::out|ios::binary));
  o << d << '\n';
  o.pop();
  o.pop();

  //now, read via boost
  filtering_istream in;
  in.push(gzip_decompressor());
  in.push(file_source("test_boostgz_out.gz",ios::in|ios::binary));
  SimData d3;
  in >> d3 >> ws;
  in.pop();
  in.pop();

  //read the boost file using the zlib routines
  gzf = gzopen("test_boostgz_out.gz","r");
  SimData d4 = read_SimData_gz(gzf);
  gzclose(gzf);

  //write it in binary
  ofstream obin("test_binary_out.bin",ios::binary);
  write_SimData_binary(obin,d);
  obin.close();

  //read it
  ifstream ibin("test_binary_out.bin",ios::binary);
  SimData d5 = read_SimData_binary(ibin);
  ibin.close();
  

  cout << d << '\n' << d2 << '\n' << d3 << '\n' << d4 << '\n' << d5 << '\n';
}
