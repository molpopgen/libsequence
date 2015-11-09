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

  This sounds cool, too: c++ concurrency book
 */
namespace Sequence
{
  std::vector<double> lHaf_t( const SimData & data, const double l, const int nthreads )
  {
    std::vector<unsigned> dcounts(data.numsites());
    
    const auto counter = [](SimData::const_site_iterator __beg,
			    const SimData::const_site_iterator __end,
			    unsigned * rv ) {
      std::vector<unsigned> temp;
      for( ; __beg < __end ; ++__beg )
	{
	  temp.emplace_back( unsigned(std::count(__beg->second.begin(),__beg->second.end(),'1')) );
	}
      std::copy(temp.begin(),temp.end(),rv);
    };
    std::vector<std::thread> vt;
    auto beg = data.sbegin();
    auto end=data.send();
    unsigned snp_per_thread = unsigned(std::ceil(double(end-beg)/double(nthreads)));
    for( beg ; beg < end ; beg += snp_per_thread)
      {
	vt.emplace_back(std::thread(counter,beg,
				    (snp_per_thread < std::distance(beg,end)) ? beg + snp_per_thread : end,
				    &dcounts[beg-data.sbegin()]));
      }
    for( auto i = 0 ; i < vt.size() ; ++i ) vt[i].join();
    std::vector<double> rv(data.size(),0.);
    vt.clear();
    const auto counter2 = []( SimData::const_data_iterator beg,
			      const SimData::const_data_iterator end,
			      const double l,
			      const unsigned * dc,
			      double * rv )
      {
	std::vector<double> temp;
	for( ; beg < end ; ++beg ) {
	  double rv_local=0.;
	  auto j = std::find_if(beg->cbegin(),beg->cend(),[](const char & ch) {
	      return ch == '1';
	    });
	  while(j != beg->cend())
	    {
	      size_t d = size_t(j-beg->cbegin());
	      rv_local += std::pow(*(dc+d),l );
	      j = std::find_if(j+1,beg->cend(),
			       [](const char & ch) {
				 return ch == '1';
			       });
	    }
	  temp.emplace_back(rv_local);
	}
	std::copy(temp.begin(),temp.end(),rv);
      };
    unsigned hap_per_thread = unsigned(std::ceil(double(data.size())/double(nthreads)));
    for(PolyTable::const_data_iterator i = data.cbegin() ; i != data.cend() ; i += hap_per_thread)
      {
	vt.emplace_back(std::thread(counter2,i,
				    (hap_per_thread < std::distance(i,data.cend())) ? (i + hap_per_thread) : data.cend(),
				    l,
				    &dcounts[0],
				    &rv[i-data.cbegin()]));
      }
    for( unsigned j = 0 ; j <vt.size() ; ++j ) vt[j].join();
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
    std::vector<double> rv(data.size(),0.);
    for(unsigned i = 0 ; i < data.size() ; ++i )
      {
	auto j = std::find_if(data[i].cbegin(),data[i].cend(),[](const char & ch) {
	    return ch == '1';
	  });
	while(j!=data[i].cend())
	  {
	    size_t d = size_t(j-data[i].cbegin());
	    rv[i] += std::pow( double(dcounts[d]),l );
	    j = std::find_if(j+1,data[i].cend(),
			     [](const char & ch) {
			       return ch == '1';
			     });
	  }
      }
    return rv;
  }
}
