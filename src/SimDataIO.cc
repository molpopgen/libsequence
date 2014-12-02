#include <Sequence/SimDataIO.hpp>
#include <Sequence/IOhelp.hpp>
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <cstdint>
#include <limits> 

using namespace std;

namespace Sequence
{
  long long write_SimData_gz( gzFile & file , const SimData & d, const bool & binary )
  {
    long long ttl = 0;
    ostringstream buffer;
    if ( binary )
      {
	write_SimData_binary(buffer,d);
	//We cannot write a buffer this big, but it is hard to imagine this situation actually happening...
	if(buffer.str().size() > numeric_limits<unsigned>::max()) return -1; 
	ttl = gzwrite( file, buffer.str().c_str(), unsigned(buffer.str().size()) );
      }
    else
      {
	buffer << d << '\n';
	//We cannot write a buffer this big, but it is hard to imagine this situation actually happening...
	if(buffer.str().size() > numeric_limits<unsigned>::max()) return -1;
	ttl = gzwrite(file, buffer.str().c_str(),unsigned(buffer.str().size()));
      }
    return ttl;
  }
  
  SimData read_SimData_gz( gzFile & file, const bool & binary )
  {
    vector<double> pos;
    vector<string>data;
    if (binary)
      {
	uint32_t nsam,nsites,nsites_i,one;
	gzread(file,&nsam,sizeof(uint32_t));
	gzread(file,&nsites,sizeof(uint32_t));
	if(!nsites) return SimData();
	pos.resize(nsites);
	gzread(file,&pos[0],nsites*sizeof(double));
	data = vector<string>(nsam,string(nsites,'0'));
	for(uint32_t i = 0 ; i < nsam ; ++i)
	  {
	    gzread(file,&nsites_i,sizeof(uint32_t));
	    for( uint32_t j = 0 ; j < nsites_i ; ++j )
	      {
		gzread(file,&one,sizeof(uint32_t));
		data[i][one]='1';
	      }
	  }
      }
    else
      {
	IOhelp::gzreaduntil()(file,[](const char & ch){ return ch == '/'; });
	string temp; //this is our buffer
	IOhelp::gzread2ws(file,temp);
	if( temp != "//" ) return SimData(); //Something is amiss with the input!
	temp.clear();
	IOhelp::gzreadws(file);
	IOhelp::gzread2ws(file,temp); //segsites: 
	IOhelp::gzreadws(file);
	temp.clear();
	IOhelp::gzread2ws(file,temp); //The number of seg sites
	unsigned long S = stoul(temp);
	IOhelp::gzreadws(file);
	IOhelp::gzread2ws(file,temp); //positions: 
	IOhelp::gzreadws(file);
	temp.clear();
	vector<double> pos;
	for( decltype(S) i = 0 ; i < S ; ++i )
	  {
	    IOhelp::gzread2ws(file,temp); //positions: 
	    pos.push_back( stod(temp) );
	    IOhelp::gzreadws(file);
	    temp.clear();
	  }
	IOhelp::gzreadws(file);
	char * haplotype = new char[S+1];
	temp.resize(S);
	vector<string> data;
	gzgets(file,&haplotype[0],int(S)+1);
	data.push_back(string(haplotype));
	char ch;
	gzread(file,&ch,sizeof(char)); //read the newline
	bool reading = true;
	while( reading )
	  {
	    int rv = gzread(file,&ch,sizeof(char)); 
	    if( rv == 0 || rv == -1 || isspace(ch) ) { 
	      delete[]haplotype;return SimData(pos,data); 
	    } //we hit eof or an error or an empty line
	    gzungetc(ch,file);
	    gzgets(file,&haplotype[0],int(S)+1);
	    data.push_back(string(haplotype));
	    rv = gzread(file,&ch,sizeof(char)); 
	    if( rv == 0 || rv == -1 || ch != '\n' ) { reading = false; }
	  }
	delete [] haplotype;
      }
    return SimData(std::move(pos),std::move(data));
  }

  void write_SimData_binary( std::ostream & o, const SimData & d )
  {
    uint32_t n = uint32_t(d.size()),numsites=uint32_t(d.numsites());
    IOhelp::writeBin(o,&n,1);
    IOhelp::writeBin(o,&numsites,1);
    for( SimData::const_pos_iterator p = d.pbegin() ; p < d.pend() ; ++p )
      {
	double pos = *p;
	IOhelp::writeBin(o,&pos,1);
      }
    for( uint32_t i=0;i<d.size();++i )
      {
	uint32_t c = uint32_t(count(d[i].begin(),d[i].end(),'1'));
	IOhelp::writeBin(o,&c,1);
	for(uint32_t j=0;j<d.numsites();++j)
	  {
	    if( d[i][j]=='1')
	      {
		IOhelp::writeBin(o,&j,1);
	      }
	  }
      }
  }

  SimData read_SimData_binary( std::istream & in )
  {
    uint32_t nsam,nsites;
    IOhelp::readBin(in,&nsam,1);
    IOhelp::readBin(in,&nsites,1);
    if( ! nsites ) { return SimData(); }
    
    vector<double> pos;
    vector<string> data;
    double p;
    for (uint32_t i = 0 ; i < nsites ; ++i )
      {
	IOhelp::readBin(in,&p,1);
	pos.push_back(p);
      }

    uint32_t nsites_i,index_i;
    for(uint32_t i = 0 ; i < nsam ; ++i )
      {
	string d(nsites,'0');
	IOhelp::readBin(in,&nsites_i,1);
	for( uint32_t j = 0 ; j < nsites_i ; ++j )
	  {
	    IOhelp::readBin(in,&index_i,1);
	    d[index_i]='1';
	  }
	data.push_back(d);
      }
    return SimData(std::move(pos),std::move(data));
  }

  long int write_SimData_binary( int fd , const SimData & d )
  {
    ostringstream o;
    write_SimData_binary( o, d );
    return write( fd, o.str().c_str(), size_t(o.str().size()) );
  }

  long int write_SimData_binary( FILE * fp, const SimData & d )
  {
    return write_SimData_binary( fileno(fp), d );
  }
} //namespace Sequence
