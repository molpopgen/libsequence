#include <Sequence/SummStats/Uhaps.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/polySiteVector.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/util/nibble.hpp>
#include <memory>
#include <utility>
#include <list>
#include <algorithm>
#include <cassert>

namespace Sequence {

  int Uhaps::populate(const polySiteVector8 & __pt)
  {
    using namespace nibble;
    if(__pt.empty()) return 0;
    //Bioloigcally, this is nsam/2...
    Seq8::size_type nsam = __pt.begin()->second.second.size();
    bool isodd = (__pt.size()%2);
    for( Seq8::size_type i = 0 ; i < nsam ; ++i )
      {
	pack8::vtype a(__pt.size()/2 + isodd);
	pack8::vtype b = (!isodd || i < nsam-1) ? a : pack8::vtype();
	bool idx = 0;
	pack8::vtype::size_type j = 0;
	for( const auto & p : __pt )
	  {
	    if( p.second.second.size() != nsam )
	      {
	     	ustrings.clear();
	     	ustring_itrs.clear();
	     	return -1;
	      }
	    if ( !idx )
	    {
	      //read hi, write hi a
	      writehi(a[j],readhi(p.second.second[i]));
	      //read lo, write hi b
	      if(!b.empty())
		{
		  writehi(b[j],readlo(p.second.second[i]));
		}
	    }
	    else 
	      {
	    	//read hi, write lo a
		writelo(a[j],readhi(p.second.second[i]));
	    	//read lo, write lo b
	    	if(!b.empty())
		  {
		    writelo(b[j],readlo(p.second.second[i]));
		  }
	      }
	    if(idx) ++j;
	    idx = !idx;
	  }
	auto itr = std::find_if( ustrings.cbegin(),
				 ustrings.cend(),
				 [&a](const Seq8 & __s8) {
				   return a == __s8.second;
				 } );
	if ( itr == ustrings.cend() )
	  {
	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(), Seq8(std::move(unsigned(2*a.size()-isodd)),std::move(a))) );
	  }
	else
	  {
	    ustring_itrs.emplace_back(std::move(itr));
	  }

	if(! b.empty() )
	  {
	     itr = std::find_if( ustrings.cbegin(),
				 ustrings.cend(),
				 [&b](const Seq8 & __s8) {
				   return b == __s8.second;
				 } );
	     if ( itr == ustrings.cend() )
	       {
		 ustring_itrs.emplace_back( ustrings.insert(ustrings.end(), Seq8(std::move(unsigned(2*b.size()-isodd)),std::move(b))) );
	       }
	     else
	       {
		 ustring_itrs.emplace_back(std::move(itr));
	       }
	  }
      }
    return 0;
  }

  int Uhaps::populate(const PolyTable & __pt)
  {
    if(__pt.empty()) return 0;
    for( const auto & s : __pt )
      {
	Seq8 t(s,dna_poly_alphabet);
	auto itr = std::find(ustrings.cbegin(),ustrings.cend(),t);
	if( itr == ustrings.end() )
	  {
	    ustring_itrs.emplace_back( ustrings.insert(ustrings.end(),std::move(t)) );
	  }
	else
	  {
	    ustring_itrs.emplace_back(std::move(itr));
	  }
      }
    return 0;
  }

  Uhaps::Uhaps( const polySiteVector8 & __pt) : ustrings ( Uhaps::ustring_ctr_t() ),
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
}
