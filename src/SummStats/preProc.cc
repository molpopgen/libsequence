#include <Sequence/preProc.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/Ptable.hpp>
#include <Sequence/stateCounter.hpp>
#include <memory>
#include <utility>
#include <list>
#include <algorithm>

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
    void populate( const Ptable & ps,
		   const bool & haveAncStates, const std::string::size_type & anc,
		   std::vector<stateCounter> & states,
		   std::vector<std::pair<bool,stateCounter> > & dstates);
  };

  void populate( const Ptable & ps,
		 const bool & haveAncStates, const std::string::size_type & anc,
		 std::vector<stateCounter> & states,
		 std::vector<std::pair<bool,stateCounter> > & dstates)
  {
    std::for_each( ps.cbegin(), ps.cend(),
		   [&states,&dstates,&haveAncStates,&anc]( const polymorphicSite & __ps )
		   {
		     stateCounter a,d; //all && derived states, resp.
		     std::string::size_type i = 0;
		     bool ancValid = haveAncStates;
		     for( auto itr = __ps.second.cbegin() ; itr != __ps.second.cend() ; ++itr,++i )
		       {
			 a(*itr);
			 if( haveAncStates )
			   {
			     if( i == anc ) 
			       {
				 if (std::toupper(*itr) == 'N' || *itr == '-' )
				   ancValid = false;
			       }
			     else d(*itr);
			   }
		       }
		     states.emplace_back(std::move(a));
		     if( haveAncStates )
		       dstates.emplace_back( std::move(std::make_pair(std::move(ancValid),std::move(d))) );
		   } );
  }

  int Uhaps::populate(const Ptable & __pt)
  {
    if(__pt.empty()) return 0;
    std::string::size_type nsam = __pt.begin()->second.size();
    for( std::string::size_type i = 0 ; i < nsam ; ++i )
      {
	std::string hap;
	for( const auto & p : __pt )
	  {
	    if( p.second.size() != nsam )
	      {
		ustrings.clear();
		ustring_itrs.clear();
		return -1;
	      }
	    hap += p.second[i];
	  }
	auto itr = std::find(ustrings.cbegin(),ustrings.cend(),hap);
	if( itr == ustrings.cend() ) //new string
	  {
	    //store the pointer and add string
	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(),hap) );
	  }
	else
	  {
	    //store the pointer
	    ustring_itrs.push_back(itr);
	  }
      }
    return 0;
  }

  int Uhaps::populate(const PolyTable & __pt)
  {
    if(__pt.empty()) return 0;
    std::string::size_type len = __pt[0].size();
    for( auto s : __pt )
      {
	if( s.size() != len ) 
	  {
	    ustrings.clear();
	    ustring_itrs.clear();
	    return -1;
	  }
	auto itr = std::find(ustrings.cbegin(),ustrings.cend(),s);
	if( itr == ustrings.cend() ) //new string
	  {
	    //store the pointer and add string
	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(),s) );
	  }
	else
	  {
	    //store the pointer
	    ustring_itrs.push_back(itr);
	  }
      }
    return 0;
  }

  Uhaps::Uhaps( const Ptable & __pt) : ustrings ( Uhaps::ustring_ctr_t() ),
				       ustring_itrs (Uhaps::ustring_itr_ctr_t() )
  {
    if( populate(__pt) == -1 )
      {
	throw SeqException("Sequence::Uhaps error.  Invalid data encountered");
      }
  }

  Uhaps::Uhaps( const PolyTable & __pt) : ustrings ( Uhaps::ustring_ctr_t() ),
					  ustring_itrs (Uhaps::ustring_itr_ctr_t() )
  {
    if( populate(__pt) == -1 )
      {
	throw SeqException("Sequence::Uhaps error.  Invalid data encountered");
      }
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
  
  preProc::preProc( Ptable && __pt, 
		    const bool & haveAncStates,
		    const std::string::size_type & anc) : __impl(new preProcImpl()),
							  ptable(std::move(__pt)),
							  uhaps(__pt),
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
}
