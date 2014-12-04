#include <Sequence/SimData.hpp>
#include <algorithm>
#include <cmath>

namespace
{
  double __nSLsum(const unsigned & core,
                  const Sequence::SimData & d,
                  const std::vector<size_t> & coretype)
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
		s += unsigned(std::distance(left.first.base(),right.first)) + 1;
		++nc;
	      }
	  }
      }
    return double(s)/double(nc);
  }
}

namespace Sequence 
{
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
}
