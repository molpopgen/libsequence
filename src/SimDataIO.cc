#include <Sequence/SimDataIO.hpp>
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

namespace Sequence
{
  long long write_SimData_gz( gzFile & file , const SimData & d )
  {
    long long ttl = 0;
    int w = gzwrite( file, "//\nsegsites: ", 13 );
    if(w < 0) { return w; }
    ttl += w;

    w = gzprintf( file, "%u\n", d.numsites() );
    if ( w < 0 ) { return w; }
    ttl += w;

    if (! d.empty() )
      {
	int w = gzwrite( file, "positions: ", 11 );
	SimData::const_pos_iterator pi = d.pbegin();
	for( ; pi < d.pend()-1 ; ++pi )
	  {
	    w = gzprintf( file, "%lf ", *pi );
	    if ( w < 0 ) { return w; }
	    ttl += w;
	  }
	w = gzprintf( file, "%lf\n", *pi );
	if ( w < 0 ) { return w; }
	ttl += w;
      
	for(unsigned ind = 0 ; ind < d.size() ; ++ind)
	  {
	    w = gzprintf( file, "%s\n", d[ind].c_str() );
	    if ( w < 0 ) { return w; }
	    ttl += w;
	  }
      }
    return ttl;
  }

  SimData read_SimData_gz( gzFile & file )
  {
    char * buffer1 = new char[13];
    char * endptr;

    gzgets(file,buffer1,100);
    if( string(buffer1) != string("//\n") ) 
      {
	delete[] buffer1;
	return SimData();
      }
    gzgets(file,buffer1,10);
    string buffer2;
    int ch = gzgetc(file);
    while( char(ch) != '\n' )
      {
	if ( !isspace(ch) )
	  {
	    buffer2 += ch;
	  }
	ch = gzgetc(file);
      }

    unsigned segsites = unsigned(strtoul(buffer2.c_str(), &endptr, 10));

    if( segsites )
      {
	gzgets(file,buffer1,12);
	vector<double> pos;
	buffer2.clear();
	ch = gzgetc(file);
	while( char(ch) != '\n' )
	  {
	    if( ! isspace(ch) )
	      {
		buffer2 += char(ch);
	      }
	    else
	      {
		pos.push_back(strtod(buffer2.c_str(),&endptr));
		buffer2.clear();
	      }
	    ch = gzgetc(file);
	  }
	pos.push_back(strtod(buffer2.c_str(),&endptr));
	//read the positions
	vector<string> data;
	ch = gzgetc(file);
	char * buffer3 = new char[segsites];
	while( char(ch) != '/' && ch != EOF )
	  {
	    gzungetc(ch,file);
	    gzgets( file,&buffer3[0], segsites+2);
	    string temp(&buffer3[0]);
	    data.push_back( string(temp.begin(),temp.end()-1) );
	    ch = gzgetc(file);	
	  }
	delete [] buffer3;
	if(char(ch) == '/')
	  {
	    gzungetc(ch,file);
	  }
	delete [] buffer1;
	return SimData(pos,data);
      }
    delete [] buffer1;
    return SimData();
  }

  void write_SimData_binary( std::ostream & o, const SimData & d )
  {
    unsigned n = d.size(),numsites=d.numsites();
    o.write( reinterpret_cast<char *>(&n),sizeof(unsigned) );
    o.write( reinterpret_cast<char *>(&numsites),sizeof(unsigned) );
    for( SimData::const_pos_iterator p = d.pbegin() ; p < d.pend() ; ++p )
      {
	double pos = *p;
	o.write( reinterpret_cast<char*>(&pos),sizeof(double) );
      }
    for( unsigned i=0;i<d.size();++i )
      {
	unsigned c = count(d[i].begin(),d[i].end(),'1');
	o.write(reinterpret_cast<char*>(&c),sizeof(unsigned));
	for(unsigned j=0;j<d.numsites();++j)
	  {
	    if( d[i][j]=='1')
	      {
		o.write( reinterpret_cast<char *>(&j),sizeof(unsigned) );
	      }
	  }
      }
  }

  SimData read_SimData_binary( std::istream & in )
  {
    unsigned nsam,nsites;
    in.read ( reinterpret_cast<char *>(&nsam), sizeof(unsigned) );
    in.read ( reinterpret_cast<char *>(&nsites), sizeof(unsigned) );
    if( ! nsites ) { return SimData(); }
    
    vector<double> pos;
    vector<string> data;
    double p;
    for (unsigned i = 0 ; i < nsites ; ++i )
      {
	in.read( reinterpret_cast<char *>(&p),sizeof(double) );
	pos.push_back(p);
      }

    unsigned nsites_i,index_i;
    for(unsigned i = 0 ; i < nsam ; ++i )
      {
	string d(nsites,'0');
	in.read( reinterpret_cast<char *>(&nsites_i), sizeof(unsigned) );
	for( unsigned j = 0 ; j < nsites_i ; ++j )
	  {
	    in.read( reinterpret_cast<char *>(&index_i), sizeof(unsigned) );
	    d[index_i]='1';
	  }
	data.push_back(d);
      }
    return SimData(pos,data);
  }

  int write_SimData_binary( int fd , const SimData & d )
  {
    ostringstream o;
    write_SimData_binary( o, d );
    return write( fd, o.str().c_str(), o.str().size() );
  }

  int write_SimData_binary( FILE * fp, const SimData & d )
  {
    return write_SimData_binary( fileno(fp), d );
  }
} //namespace Sequence
