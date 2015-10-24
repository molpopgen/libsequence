#include <Sequence/SummStats/lHaf.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <thread>

/*
  TODO:  A thread_site and thread_hap wrapper that apply some fxn to a site
  or a hap, using <thread>, filling up a vector of return values.

  Need to read/remind myself of this:

  http://www.aristeia.com/TalkNotes/ACCU2011_CPUCaches.pdf

  This sounds cool, too:
 */
namespace Sequence
{
  std::vector<double> lHaf_t( const SimData & data, const double l, const int nthreads )
  {
    std::vector<unsigned> dcounts(data.numsites());
    
    //std::function<void(const SimData::const_site_iterator& , unsigned &)>
    const auto counter = [](const SimData::const_site_iterator& __s , unsigned & i ) {
      i=unsigned(std::count(__s->second.begin(),__s->second.end(),'1'));
    };
    
    for( auto i = data.sbegin() ; i != data.send() ; )
      {
	std::vector<std::thread> vt;
	for(int j = 0 ; i != data.send() && j < nthreads ; ++i,++j)
	  {
	    vt.emplace_back(std::thread(std::cref(counter),
					i,std::ref(dcounts[unsigned(i-data.sbegin())])));
	  }
	for( unsigned j = 0 ; j < unsigned(nthreads) ; ++j ) vt[j].join();
      }
    std::vector<double> rv(data.size());
    for(auto i = data.cbegin() ; i != data.cend() ; )
      {
	std::vector<std::thread> vt;
	for(int j = 0 ; i != data.end() && j < nthreads ; ++i,++j )
	  {
	    vt.emplace_back(std::thread([&rv,&dcounts,&j,&data,l](const SimData::const_data_iterator & __hap) {
		  double rv_local = 0.;
		  auto j = std::find_if(__hap->cbegin(),__hap->cend(),[](const char & ch) {
		      return ch == '1';
		    });
		while(j != __hap->cend())
		  {
		    size_t d = size_t(j-__hap->cbegin());
		    rv_local += std::pow( double(dcounts[d]),l );
		    j = std::find_if(j+1,__hap->cend(),
				     [](const char & ch) {
				       return ch == '1';
				     });
		  }
		rv[unsigned(__hap-data.begin())] = rv_local;
		},i));
	  } 
      }
    return rv;
  }

  
  std::vector<double> lHaf( const SimData & data, const double l )
  {
    //Get derived mutation frequency counts per site
    std::vector<unsigned> dcounts;
    std::for_each( data.sbegin(), data.send(),
		   [&dcounts]( const polymorphicSite & p ) {
		     dcounts.push_back( unsigned(std::count(p.second.begin(),p.second.end(),'1')) );
		   } );
    //Get the values for each element in the data
    std::vector<double> rv(data.numsites(),0.);
    for(unsigned i = 0 ; i < data.size() ; ++i )
      {
	auto j = std::find_if(data[i].cbegin(),data[i].cend(),[](const char & ch) {
	    return ch == '1';
	  });
	while(j!=data[i].cend())
	  {
	    size_t d = size_t(j-data[i].cbegin());
	    rv[d] += std::pow( double(dcounts[d]),l );
	    j = std::find_if(j+1,data[i].cend(),
			     [](const char & ch) {
			       return ch == '1';
			     });
	  }
      }
    return rv;
  }
}
