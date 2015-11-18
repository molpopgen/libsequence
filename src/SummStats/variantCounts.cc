#include <Sequence/SummStats/variantCounts.hpp>
namespace Sequence
{
  template<>
  variantCounts<SimData> make_variantCounts(const SimData & sd, bool , unsigned, char  )
  {
    return variantCounts<SimData>(sd,false);
  }
 }
