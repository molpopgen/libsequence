// Code for the -*- C++ -*- namespace Sequence
#include <type_traits>
#include <algorithm>
#include <Sequence/stateCounter.hpp>
namespace Sequence
{
  template<typename T> T removeGaps( const T & t, const bool skipAnc, const unsigned anc,const char gapchar)
  {
    auto remover = [](const stateCounter & sc) { return !sc.gap; };
    return removeColumns(t,remover,skipAnc,anc,gapchar);
  }

  template<typename T,
	   typename F> T removeColumns(const T & t,
				       const F & f,
				       const bool skipAnc,
				       const unsigned anc,
				       const char gapchar)
  {
    static_assert( std::is_base_of<PolyTable,T>::value,
		   "T must be derived from Sequence::PolyTable" );
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
	    sc = std::for_each(std::begin(__c.second)+anc+1,std::end(__c.second),sc);
	    if (f(sc)) columns.push_back(__c);
	  });
      }
    return T(columns.cbegin(),columns.cend());
  }
  
  template<typename T> T removeInvariantPos(const T & t, const bool skipAnc,
					    const unsigned anc,
					    const char gapchar)
  {
    auto remover = [](const stateCounter & sc) {return sc.nStates()>1;};
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

  template<typename T> T removeMultiHits(const T & t, const bool skipAnc,
					 const unsigned anc,
					 const char gapchar)
  {
    auto remover = [](const stateCounter & sc) { return sc.nStates()<=2; };
    return removeColumns(t,remover,skipAnc,anc,gapchar);
  }

  template<typename T> T polyTableToBinary(const T & t, const unsigned ref, const char gapchar)
  {
    T t2 = removeMultiHits(t,true,ref,gapchar);
    if (t2.empty()) return t2;
    std::vector<std::string> bd(t2.size(),std::string(t2.numsites(),'0'));
    auto updater = [&t2,&bd,&ref,gapchar](unsigned f, unsigned l)
      {
	for( ; f < l ; ++f )
	  {
	    for(unsigned i = 0 ; i < t2[f].size() ; ++i)
	      {
		if(t2[f][i]!=t2[ref][i]) bd[f][i]='1';
	      }
	  }
      };
    updater(0,ref);
    updater(ref+1,t2.size());
    return T(std::vector<double>(t2.pbegin(),t2.pend()),
	     std::move(bd));
  }
}
