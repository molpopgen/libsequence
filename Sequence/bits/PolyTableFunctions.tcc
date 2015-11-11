// Code for the -*- C++ -*- namespace Sequence
#include <type_traits>
#include <algorithm>
#include <Sequence/stateCounter.hpp>
namespace Sequence
{
  template<typename T> T removeGaps( const T & t, const char gapchar)
  {
    std::vector<typename T::column_t> columns;
    std::for_each( t.sbegin(),t.send(),
		   [&columns,gapchar](const typename T::column_t & __c) {
		     if (__c.second.find(gapchar)==std::string::npos) {
		       columns.push_back(__c);
		     }
		   });
    return T(columns.cbegin(),columns.cend());
  }

  
  template<typename T> T removeInvariantPos(const T & t, const bool skipAnc,
					    const unsigned anc,
					    const char gapchar)
  {
    std::vector<typename T::column_t> columns;
    if(!skipAnc)
      {
	std::for_each( t.sbegin(),t.send(),[&columns,gapchar](const typename T::column_t & __c) {
	    auto sc = std::for_each(std::begin(__c.second),std::end(__c.second),stateCounter(gapchar));
	    if (sc.nStates()==1 && !sc.ndna) columns.push_back(__c);
	  });
      }
    else
      {
	std::for_each( t.sbegin(),t.send(),[&columns,anc,gapchar](const typename T::column_t & __c) {
	    auto sc = std::for_each(std::begin(__c.second),std::begin(__c.second)+anc,stateCounter(gapchar));
	    std::for_each(std::begin(__c.second)+anc+1,std::end(__c.second),stateCounter(gapchar));
	    if (sc.nStates()==1 && !sc.ndna) columns.push_back(__c);
	  });
      }
    return T(columns.cbegin(),columns.cend());
  }
}
