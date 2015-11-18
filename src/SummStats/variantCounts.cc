#include <Sequence/SummStats/variantCounts.hpp>
namespace Sequence
{
  variableSiteData::variableSiteData( double pos_, bool haveAnc,
				      stateCounter counts_,stateCounter dcounts_) :
    pos(pos_),
    counts(std::move(counts_)),
    dcounts(std::move(dcounts_)),
    ancState(haveAnc)
  {
  }
  
  template<>
  variantCounts<SimData> make_variantCounts(const SimData & sd, bool , unsigned, char  )
  {
    return variantCounts<SimData>(sd,false);
  }
 }
