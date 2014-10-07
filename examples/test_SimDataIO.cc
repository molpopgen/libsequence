#include <Sequence/SimDataIO.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Sequence;

void print_problems( const SimData & d,
		     const SimData & d2 );

int main( int argc, char ** argv )
{
  SimData d;
  while(!cin.eof())
    {
      cin >> d >> ws;

      cerr << "Writing in gzip\n";
      gzFile gzf = gzopen("test_zlib_out.gz","w");
      write_SimData_gz(gzf, d);
      gzclose(gzf);

      cerr << "Reading from gzip\n";
      //now, try to read it
      gzf = gzopen("test_zlib_out.gz","r");
      SimData d2 = read_SimData_gz(gzf);
      gzclose(gzf);

      cerr << "Writing in binary\n";
      //write it in binary
      ofstream obin("test_binary_out.bin",ios::binary);
      write_SimData_binary(obin,d);
      obin.close();

      cerr << "Reading from binary\n";
      //read it
      ifstream ibin("test_binary_out.bin",ios::binary);
      SimData d3 = read_SimData_binary(ibin);
      ibin.close();
  
      cerr << "Writing to binary + gzip\n";
      gzf = gzopen("test_zlib_out.bin.gz","bw");
      write_SimData_gz(gzf, d, true);
      gzclose(gzf);

      cerr << "Reading from binary + gzip\n";
      //now, try to read it
      gzf = gzopen("test_zlib_out.bin.gz","r");
      SimData d4 = read_SimData_gz(gzf,true);
      gzclose(gzf);

      if( d != d2 )
	{
	  cerr << "Error: d != d2\n";
	  print_problems(d,d2);
	}
      if(d != d3)
	{
	  cerr << "Error: d != d3\n";
	  print_problems(d,d3);
	}
      if(d != d4)
	{
	  cerr << "Error: d != d4\n";
	  print_problems(d,d4);
	}
    }
}

void print_problems( const SimData & d,
		     const SimData & d2 )
{
  for( unsigned i = 0 ; i < d.numsites() ; ++i )
    {
      if( d.position(i) != d2.position(i) )
	{
	  cerr << "Position " << i << ": " << d.position(i) << ' ' << d2.position(i) << '\n';
	}
      for( unsigned i = 0 ; i < d.size() ; ++i )
	{
	  if( d[i] != d2[i] )
	    {
	      cerr << "Haplotype " << i << ": " << d[i] << "\n\n" << d2[i] << '\n';
	    }
	}
    }
}
