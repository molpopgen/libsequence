#include <algorithm>
#include <set>
#include <string>
#include <cmath>
#include <Sequence/SummStats/Garud.hpp>

using namespace std;

namespace Sequence
{
  GarudStats::GarudStats() : H1(1.),
			     H12(std::numeric_limits<double>::quiet_NaN()),
			     H2H1(std::numeric_limits<double>::quiet_NaN())
  {
  }
  
  GarudStats::GarudStats(const double __h1,
			 const double __h12,
			 const double __h2h1) : H1(__h1),
						H12(__h12),
						H2H1(__h2h1)
  {
  }
  
  GarudStats H1H12(const SimData & d)
  {
    if( d.empty() ) return GarudStats();
    set<string> uhaps(d.begin(),d.end());
    vector<double> hapcounts;
    //GarudStats G;
    double H1 = 0.,H12=0.,H2H1=0.;
    
    for_each(uhaps.cbegin(),uhaps.cend(),
	     [&](const string & hap) {
	       unsigned hapcount = unsigned(count(d.begin(),d.end(),hap));
	       //Unbiased calc. of homozygosity in finite sample
	       H1 += double(hapcount)*double(hapcount-1)/(double(d.size())*double(d.size()-1));
	       hapcounts.push_back(double(hapcount));
	     });
    sort(hapcounts.begin(),hapcounts.end(),
	 std::bind(greater<double>(),std::placeholders::_1,std::placeholders::_2));
    H12 = H1 + 2.*hapcounts[0]*hapcounts[1]/std::pow(double(d.size()),2.);
    H2H1 = (H1-double(hapcounts[0]*(hapcounts[0]-1))/double(d.size()*(d.size()-1)))/H1;
    return GarudStats(H1,H1,H2H1);
  }
}
