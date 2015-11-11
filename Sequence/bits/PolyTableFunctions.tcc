// Code for the -*- C++ -*- namespace Sequence
#include <type_traits>
#include <algorithm>
#include <Sequence/stateCounter.hpp>
namespace Sequence
{
  template<typename T> T removeGaps( const T & t, const bool skipAnc, const unsigned anc,const char gapchar)
  {
    std::vector<typename T::column_t> columns;
    std::for_each( t.sbegin(),t.send(),
		   [&columns,skipAnc,anc,gapchar](const typename T::column_t & __c) {
		     auto p = __c.second.find(gapchar);
		     if (p == std::string::npos || (skipAnc && unsigned(p) == anc))
		       {
			 columns.push_back(__c);
		       }
		   }
		   );
    return T(columns.cbegin(),columns.cend());
  }

  template<typename T,
	   typename F> T removeColumns(const T & t,
				       const F & f,
				       const bool skipAnc,
				       const unsigned anc,
				       const char gapchar)
  {
    std::vector<typename T::column_t> columns;
    if(!skipAnc)
      {
	std::for_each( t.sbegin(),t.send(),[&columns,&f,gapchar](const typename T::column_t & __c) {
	    auto sc = std::for_each(std::begin(__c.second),std::end(__c.second),stateCounter(gapchar));
	    if (f(sc)) columns.push_back(__c);
	  });
      }
    else
      {
	std::for_each( t.sbegin(),t.send(),[&columns,&f,anc,gapchar](const typename T::column_t & __c) {
	    auto sc = std::for_each(std::begin(__c.second),std::begin(__c.second)+anc,stateCounter(gapchar));
	    std::for_each(std::begin(__c.second)+anc+1,std::end(__c.second),stateCounter(gapchar));
	    if (f(sc)) columns.push_back(__c);
	  });
      }
    return T(columns.cbegin(),columns.cend());
  }
  
  template<typename T> T removeInvariantPos(const T & t, const bool skipAnc,
					    const unsigned anc,
					    const char gapchar)
  {
    auto remover = [](const stateCounter & sc) { return sc.nStates()==1 && !sc.ndna; };
    return removeColumns(t,remover,skipAnc,anc,gapchar);
  }

  template<typename T> T removeAmbiguous(const T & t, const bool skipAnc,
					 const unsigned anc,
					 const char gapchar)
  {
    auto remover = [](const stateCounter & sc) { return !sc.ndna; };
    return removeColumns(t,remover,skipAnc,anc,gapchar);
  }

  template<typename T> T removeMissing(const T & t, const bool skipAnc,
				       const unsigned anc,
				       const char gapchar)
  {
    auto remover = [](const stateCounter & sc) { return !sc.n; };
    return removeColumns(t,remover,skipAnc,anc,gapchar);
  }
}
