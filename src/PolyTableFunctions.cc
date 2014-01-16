/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/SeqProperties.hpp>
#include <Sequence/PolyTable.hpp>
#include <algorithm>
#include <set>
#include <cctype>

namespace Sequence
{
  bool containsCharacter(const PolyTable * t,
			 const char & ch)
  /*!
    \return true if \a t contains \a ch, false otherwise
  */
  {
    for( PolyTable::const_data_iterator itr = t->begin() ;
	 itr < t->end() ;
	 ++itr )
      {
	if ( itr->find(ch) != std::string::npos )
	  {
	    return true;
	  }
      }
    return false;
  }

  void fillIn(PolyTable * t,
	      const unsigned & refseq,
	      const char & identical)
  /*!
    Sometimes polymorphism data contain a special character
    that means that a particular state is identical to a
    reference sequence in the data.  This function replaces
    that character with the state of the reference sequence.
    \param t a PolyTable
    \param refseq the index of the reference sequence
    \param identical the character used to represent identity to the refseq
   */
  {
    if(t->empty()) return;
    for( PolyTable::data_iterator itr = t->begin() ;
	 itr < t->end() ;
	 ++itr )
      {
	std::string::iterator si1=itr->begin(),si2;
	while( (si2 = std::find(si1,itr->end(),identical)) != itr->end() )
	  {
	    size_t offset = si2 - itr->begin();
	    *si2 = *((t->begin()+refseq)->begin()+offset);
	    si1 = si2+1;
	  }
      }
  }

  void addIdentityChar(PolyTable *t,
		       const unsigned & refseq,
		       const char & identical)
  /*!
    Fill in a PolyTable with characters representing identity to some
    reference sequence ("refseq") in the data.
  */
  {
    if(t->empty()) return;
    const PolyTable::const_data_iterator beg = t->begin();
    for( PolyTable::data_iterator itr = t->begin() ;
	 itr < t->end() ;
	 ++itr)
      {
	size_t offset_data = itr-t->begin();
	if ( offset_data != refseq )
	  {
	    for (unsigned i = 0 ; i < (*itr).length() ; ++i)
	      {
		if ( (*itr)[i] == (*(beg+refseq))[i] )
		  {
		    (*itr)[i] = identical;
		  }
	      }
	  }
      }
  }

  void 
  RemoveGaps(PolyTable * table,
	     const char & gapchar) 

  /*!
    Removes all positions containing \a gapchar from \a table
    \

  */
  {
    if(table->empty()) return;
    std::vector<double> newpos;
    std::vector<std::string> newdata(table->size());
    for (unsigned site=0;site<table->numsites();++site)
      {
	bool gapped = false;
	for(unsigned ind=0 ; ind < table->size() ; ++ind)
	  {
	    if ((*table)[ind][site] == gapchar)
	      {
		gapped = true;
		ind = unsigned(table->size());
	      }
	  }
	if (gapped == false)
	  {
	    newpos.push_back(table->position(site));
	    for(unsigned ind=0 ; ind < table->size() ; ++ind)
	      {
		newdata[ind] += (*table)[ind][site];
	      }
	  }
      }
    if (table->assign(&newpos[0],newpos.size(),&newdata[0],newdata.size()) == false)
      {
	throw (SeqException("Sequence::RemoveGaps -- error: could not assign data to object"));
      }
  }

  void RemoveInvariantColumns(PolyTable * table,
			      const bool & skipOutgroup,
			      const unsigned & outgroup) 

  /*!
    Goes through the data and removes any columns that
    contain only 1 state.  There is an option to ignore
    the state of the outgroup in this calculation.
    \

  */
  {
    if(table->empty()) return;
    std::vector<double> _newpos;
    std::vector<std::string> _newdata(table->size());

    for(unsigned site=0;site< table->numsites();++site)
      {
	std::set<char> sc;
	for(unsigned ind=0;ind<table->size();++ind)
	  {
	    if ( !skipOutgroup || (skipOutgroup && ind != outgroup) )
	      {
		char uch =  char(std::toupper((*table)[ind][site]));
		if ( uch != 'N' )
		  sc.insert(uch);
	      }
	  }
	if(sc.size() > 1)
	  {
	    _newpos.push_back(table->position(site));
	    for(unsigned ind=0;ind<table->size();++ind)
	      {
		_newdata[ind] += (*table)[ind][site];
	      }
	  }
      }
    if (table->assign(&_newpos[0],_newpos.size(),&_newdata[0],_newdata.size()) == false)
      {
	throw (SeqException("Sequence::RemoveGaps -- error: could not assign data to object"));
      }
  }

  bool PolyTableValid(const PolyTable * table)
  /*!
    \return true if the following conditions are met : First, the length of every
    row in the table is equal to the length of the vector of positions.  Second,
    all of the characters in the data are members of the set {A,G,C,T,N,-} 
    (case-insensitive).

    This function is useful if you play around with PolyTable objects in
    non-const contexts, or read them in from files and need to check that
    the data are compatible with other routines in this library.  This 
    routine can be thought of as a PolyTable equivalent to Alignment::validForPolyAnalysis,
    which works on ranges of Sequence::Seq objects.
   */
  {
    for ( PolyTable::const_data_iterator itr = table->begin() ;
	  itr < table->end() ; 
	  ++itr )
      {
	if ( (std::find_if(itr->begin(),itr->end(),invalidPolyChar()) != itr->end())
	     || ( itr->length() != table->numsites() ) )
	  {
	    return false;
	  }
      }
    return true;
  }
}
