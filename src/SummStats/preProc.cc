#include <Sequence/preProc.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/Ptable.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/util/nibble.hpp>
#include <memory>
#include <utility>
#include <list>
#include <algorithm>
#include <iostream>
#include <cassert>
/*
  Idea:
  3 classes:
  1. statePreproc
  2. derstatePreproc
  3. HaploPreproc

  The main class contains one of each of them.
  Probably mutable.

  Simplifies iteration down to beg/end, etc.

  template base class to enforce some things?

  Maybe too complex?

  Take a cue from the STL and have "tags".  This forces all Summ stats to be templates, though...
*/
namespace Sequence {

  struct preProcImpl
  {
    void populate( const Ptable8 & ps,
		   const bool & haveAncStates, const std::string::size_type & anc,
		   std::vector<stateCounter> & states,
		   std::vector<std::pair<bool,stateCounter> > & dstates);
  };

  void preProcImpl::populate( const Ptable8 & ps,
			      const bool & haveAncStates, const std::string::size_type & anc,
			      std::vector<stateCounter> & states,
			      std::vector<std::pair<bool,stateCounter> > & dstates)
  {
    using namespace nibble;
    std::for_each( ps.cbegin(), ps.cend(),
    		   [&states,&dstates,&haveAncStates,&anc]( const poly8site & __ps )
    		   {
    		     stateCounter a,d; //all && derived states, resp.
		     std::for_each( __ps.second.begin(),
				    __ps.second.end(),
				    [&a]( const poly8::itype & __i )
				    {
				      if( readhi(__i) >= NOTPOLYCHAR ) abort;
				      if( readlo(__i) >= NOTPOLYCHAR ) abort;
				      a( dna_poly_alphabet[ readhi(__i) ] );
				      a( dna_poly_alphabet[ readlo(__i) ] );
				    } );
		     states.emplace_back(std::move(a));
    		   } );
  }

  int Uhaps::populate(const Ptable8 & __pt)
  {
    using namespace nibble;
    if(__pt.empty()) return 0;
    poly8::vtype::size_type nsam = __pt.begin()->second.size();
    std::cerr << "nsam = " << nsam << '\n';
    auto X = POLYEOS;
    poly8::itype EMPTY = 0;
    writehi(EMPTY,poly8::itype(POLYEOS));
    //std::cerr << "EMPTY = " << int(EMPTY) << '\n';
    writelo(EMPTY,poly8::itype(POLYEOS));
    //std::cerr << "EMPTY = " << int(EMPTY) << '\n';
     for( poly8::vtype::size_type i = 0 ; i < nsam ; ++i )
       {
	 poly8::vtype a(__pt.size()/2+1,EMPTY),b(__pt.size()/2+1,EMPTY);
	 bool idx = 1;
	 poly8::vtype::size_type j = 0;
	 for( const auto & p : __pt )
	   {
	     // std::cerr << poly8::vtype2dna(p.second) << ' ' << p.second.size() << '\n';
	     if( p.second.size() != nsam )
	       {
		 ustrings.clear();
		 ustring_itrs.clear();
		 return -1;
	       }
	     if ( idx )
	       {
		 // std::cerr << "!idx: "<< int(readhi(p.second[i])) << ' ' << int(readlo(p.second[i])) << '\n';
		 //read hi, write hi a
		 writehi(a[j],readhi(p.second[i]));
		 //read lo, write hi b
		 writehi(b[j],readlo(p.second[i]));
	       }
	     else 
	       {
		 //std::cerr << "idx: "<< int(readhi(p.second[i])) << ' ' << int(readlo(p.second[i])) << '\n';
		 //read hi, write lo a
		 writelo(a[j],readhi(p.second[i]));
		 //read lo, write lo b
		 writelo(b[j],readlo(p.second[i]));
	       }
	     //std::cerr << i << ' ' << j << '\n';
	     idx = !idx;
	     if(idx) ++j;
	   }
	 // std::cerr << "a = ";
	 //for( auto __a : a ) std::cout << int(__a) << ' ';
	 //std::cerr << ' ' << poly8::vtype2dna(std::cref(a)) << '\n';
	 if( !(readhi(a[0]) == poly8::itype(X)) )
	   {
	     auto itr = std::find(ustrings.cbegin(),ustrings.cend(),a);
	     if( itr == ustrings.cend() ) //new hap
	       {
		 ustring_itrs.emplace_back( ustrings.insert(ustrings.end(), std::move(a)) );
	       }
	 else ustring_itrs.push_back(itr);
	   }
	 // auto bblank = std::count_if(b.begin(),b.end(),
	 // 			     [&X](const poly8::itype & i) {
	 // 			       return readhi(i) == poly8::itype(X) || readlo(i) == poly8::itype(X);
	 // 			     });
	 // if(bblank < b.size() )
	 if( !(readhi(b[0]) == poly8::itype(X) ) )
	 {
	   auto itr = std::find(ustrings.cbegin(),ustrings.cend(),b);
	   if( itr == ustrings.cend() ) //new hap
	     {
	       ustring_itrs.emplace_back( ustrings.insert(ustrings.end(), std::move(b)) );
	     }
	   else ustring_itrs.push_back(itr);
	 }
       }
     std::cerr << "Uhaps::Ptable8: "<< ustrings.size() << '\n';
    //   {
    // 	std::string hap(__pt.size(),' ');
    // 	unsigned j = 0;
    // 	for( const auto & p : __pt )
    // 	  {
    // 	    if( p.second.size() != nsam )
    // 	      {
    // 		ustrings.clear();
    // 		ustring_itrs.clear();
    // 		return -1;
    // 	      }
    // 	    hap[j++]=p.second[i];
    // 	  }
    // 	auto itr = std::find(ustrings.cbegin(),ustrings.cend(),hap);
    // 	if( itr == ustrings.cend() ) //new string
    // 	  {
    // 	    //store the pointer and add string
    // 	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(),std::move(hap)) );
    // 	  }
    // 	else
    // 	  {
    // 	    //store the pointer
    // 	    ustring_itrs.push_back(itr);
    // 	  }
    //   }
    // std::cerr << ustrings.size() << '\n';
    return 0;
  }

  int Uhaps::populate(const PolyTable & __pt)
  {
    if(__pt.empty()) return 0;
    std::string::size_type len = __pt[0].size();
    for( const auto  & s : __pt )
      {
    	if( s.size() != len ) 
    	  {
    	    ustrings.clear();
    	    ustring_itrs.clear();
    	    return -1;
    	  }
	poly8::vtype temp = poly8::dna2vtype( s );
    	auto itr = std::find(ustrings.cbegin(),ustrings.cend(),temp);
    	if( itr == ustrings.cend() ) //new string
    	  {
    	    //store the pointer and add string
    	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(),std::move(temp)) );
    	  }
    	else
    	  {
    	    //store the pointer
    	    ustring_itrs.push_back(itr);
    	  }
      }
    return 0;
  }

  int Uhaps::populate(PolyTable && __pt)
  {
    std::cerr << "move version\n";
    if(__pt.empty()) return 0;
    std::string::size_type len = __pt[0].size();
    for( auto & s : __pt )
      {
    	if( s.size() != len ) 
    	  {
    	    ustrings.clear();
    	    ustring_itrs.clear();
    	    return -1;
    	  }
	poly8::vtype temp = poly8::dna2vtype( std::move(s) );
	std::cerr << s.size() << '\n';
    	//auto itr = std::find(ustrings.cbegin(),ustrings.cend(),s);
	auto itr = std::find(ustrings.cbegin(),ustrings.cend(),temp);
    	if( itr == ustrings.cend() ) //new string
    	  {
    	    //store the pointer and add string
    	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(),std::move(temp)) );
    	  }
    	else
    	  {
    	    //store the pointer
    	    ustring_itrs.push_back(itr);
    	  }
      }
    //leave __pt in a sensible state for the calling environment...
    __pt.first.clear();
    __pt.second.clear();
    return 0;
  }

  Uhaps::Uhaps( const Ptable8 & __pt) : ustrings ( Uhaps::ustring_ctr_t() ),
					ustring_itrs (Uhaps::ustring_itr_ctr_t() )
  {
    if( populate(__pt) == -1 )
      {
	throw SeqException("Sequence::Uhaps error.  Invalid data encountered");
      }
    std::cerr << ustrings.size() << ' ' << ustring_itrs.size() << '\n';
  }

  Uhaps::Uhaps( const PolyTable & __pt) : ustrings ( Uhaps::ustring_ctr_t() ),
					  ustring_itrs (Uhaps::ustring_itr_ctr_t() )
  {
    if( populate(__pt) == -1 )
      {
	throw SeqException("Sequence::Uhaps error.  Invalid data encountered");
      }
    std::cerr << ustrings.size() << ' ' << ustring_itrs.size() << '\n';
  }

  Uhaps::Uhaps(  PolyTable && __pt) : ustrings ( Uhaps::ustring_ctr_t() ),
				      ustring_itrs (Uhaps::ustring_itr_ctr_t() )
  {
    if( populate(std::move(__pt)) == -1 )
      {
	throw SeqException("Sequence::Uhaps error.  Invalid data encountered");
      }
    std::cerr << ustrings.size() << ' ' << ustring_itrs.size() << '\n';
  }

  Uhaps::iterator Uhaps::begin() {
    return ustring_itrs.begin();
  }

  Uhaps::iterator Uhaps::end() {
    return ustring_itrs.end();
  }
  Uhaps::const_iterator Uhaps::begin() const {
    return ustring_itrs.begin();
  }

  Uhaps::const_iterator Uhaps::end() const {
    return ustring_itrs.end();
  }

  Uhaps::const_iterator Uhaps::cbegin() const {
    return ustring_itrs.cbegin();
  }

  Uhaps::const_iterator Uhaps::cend() const {
    return ustring_itrs.cend();
  }

  Uhaps::const_reference Uhaps::operator[]( const Uhaps::size_type & i ) const {
    return *ustring_itrs[i];
  }

  Uhaps::size_type Uhaps::size() const
  {
    return ustring_itrs.size();
  }

  bool Uhaps::empty() const
  {
    return ustring_itrs.empty();
  }

  /*
  preProc::preProc() : __impl(new preProcImpl()),
		       ptable(Ptable()),
		       uhaps()
  {
  }
  */
  preProc::preProc( const Ptable & __pt, 
		    const bool & haveAncStates , 
		    const std::string::size_type&  anc ) : __impl(new preProcImpl()),
							   ptable(__pt),
							   uhaps(__pt),
							   states( state_t() ),
							   dstates( dstate_t() )
  {
    __impl->populate(ptable,haveAncStates,anc,
		     states,dstates);
  }
  
  preProc::preProc( Ptable & __pt, 
		    const bool & haveAncStates,
		    const std::string::size_type & anc) : __impl(new preProcImpl()),
							  ptable(__pt),
							  uhaps(ptable),
							  states( state_t() ),
							  dstates( dstate_t() )
  {
    __impl->populate(ptable,haveAncStates,anc,
		     states,dstates);
  }
  
  preProc::preProc( const PolyTable & __pt, 
		    const bool & haveAncStates,
		    const PolyTable::size_type & anc ) : __impl(new preProcImpl()),
							 ptable(__pt),
							 uhaps(__pt),
							 states( state_t() ),
							 dstates( dstate_t() )
  {
    __impl->populate(ptable,haveAncStates,anc,
		     states,dstates);
  }

  preProc::preProc( PolyTable && __pt, 
		    const bool & haveAncStates,
		    const PolyTable::size_type & anc ) : __impl(new preProcImpl()),
							 ptable(__pt),
							 uhaps(std::move(__pt)),
							 states( state_t() ),
							 dstates( dstate_t() )
  {
    __impl->populate(ptable,haveAncStates,anc,
		     states,dstates);
  }

  preProc::~preProc(void)
  {
  }
}
