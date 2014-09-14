#include <Sequence/Ptable.hpp>
#include <Sequence/SimData.hpp>
#include <algorithm>
#include <cctype>
using namespace std;

namespace { //prototypes for non-exported functions
  //Implementation details for nSL statistic
  double __nSLsum(const unsigned & core,
                  const Sequence::SimData & d,
                  const std::vector<size_t> coretype);
}
namespace Sequence
{
  double Dij(const polymorphicSite & p, const std::vector< unsigned > & config, const unsigned & i, const unsigned & j)
  {
    unsigned rv = 0;
    unsigned N = 0;
    unsigned start1 = accumulate(config.begin(),config.begin()+i,0u),
      start2 = accumulate(config.begin(),config.begin()+j,0u);
    for( unsigned x = start1 ; x < start1 + config[i] ; ++x )
      {
	for(unsigned y = start2 ; y < start2 + config[j] ; ++y)
	  {
	    char ch1 = char(std::toupper(p.second[x])),ch2=char(std::toupper(p.second[y]));
	    if(ch1 != 'N' && ch2 != 'N')
	      {
		rv += (ch1 != ch2) ? 1u : 0u;
	      }
	    else
	      {
		++N;
	      }
	  }
      }
    return double(rv)/(double(config[i]+config[j]-N));
  }
  
  double Gmin(const Ptable & pt, const std::vector< unsigned > & config)
  {
    unsigned mdxy = numeric_limits<unsigned>::max(),sumdxy=0;
  }

  double nSL(const unsigned & core,
	     const SimData & d)
  {
    std::vector<size_t> der,anc;
    for(unsigned i=0;i<d.size();++i)
      {
	if( d[i][core] == '1' ) der.push_back(i);
	else anc.push_back(i);
      }
    //double SLDK = 2.*double(__nlSsum(core,d,der))/double(der.size()*(der.size()-1));
    //double SLAK = 2.*double(__nlSsum(core,d,anc))/double(anc.size()*(anc.size()-1));
    //return log(SLAK) - log(SLDK);
    return log(__nSLsum(core,d,anc)) - log(__nSLsum(core,d,der));
  }
}//ns Sequence

namespace{
  double __nSLsum(const unsigned & core,
                  const Sequence::SimData & d,
                  const std::vector<size_t> coretype)
  {
    unsigned s = 0,nc=0;
    for( unsigned i = 0 ; i < coretype.size() ; ++i )
      {
	for( unsigned j = i+1 ; j < coretype.size() ; ++j )
	  {
	    auto right = std::mismatch(d[coretype[i]].cbegin()+core,d[coretype[i]].cend(),
				  d[coretype[j]].cbegin()+core);
	    std::string::const_reverse_iterator ri1(d[coretype[i]].cbegin()+core),
	      ri2(d[coretype[j]].cbegin()+core);
	    auto left = std::mismatch(ri1,d[coretype[i]].crend(),ri2);
	    if(left.first != d[coretype[i]].rend() && right.first != d[coretype[i]].end())
	      {
		s += std::distance(left.first.base(),right.first) + 1;
		++nc;
	      }
	  }
      }
    return double(s)/double(nc);
  }
}
