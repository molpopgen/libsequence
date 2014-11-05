#include <Sequence/samfunctions.hpp>
#include <limits>
#include <cstdlib>
#include <algorithm>

using namespace std;

namespace {
  vector<pair<char,unsigned> > parse_cigar(const string & cigar)
  {
    vector<pair<char,unsigned> > cigar_data;
    
    string::const_iterator pibeg = find_if( cigar.begin(),cigar.end(),::isdigit);
    string::const_iterator piend = find_if( pibeg+1,cigar.end(),::isalpha );
    char * endptr;
    while( pibeg != cigar.end() )
      {
	cigar_data.push_back( make_pair( *piend, 
					 strtoul( string(pibeg,piend).c_str(),&endptr,10 ) ) );
	pibeg = find_if( piend+1, cigar.end(), ::isdigit );
	piend = find_if( pibeg+1,cigar.end(),::isalpha);
      }
    return cigar_data;
  }
}

namespace Sequence
{
  unsigned alignment_length(const samrecord & b)
  /*!
    \param b A Sequence::samrecord
    \return The sum of all M,I,D, and N elements of a cigar string
  */
  {
    unsigned sum=0;
    for( samrecord::cigar_iterator i = b.cig_begin() ;
	 i < b.cig_end() ; ++i )
      {
	switch( i->first )
	  {
	  case 'M':
	    sum += i->second;
	    break;
	  case 'I':
	    sum += i->second;
	    break;
	  case 'D':
	    sum += i->second;
	    break;
	  case 'N':
	    sum += i->second;
	    break;
	  }
      }
    return sum;
  }

  unsigned alignment_length(const bamrecord & b)
  /*!
    \param b A Sequence::bamrecord
    \return The sum of all M,I,D, and N elements of a cigar string
  */
  {
    auto cigdata = parse_cigar(b.cigar());
    unsigned sum = 0;
    std::for_each(cigdata.cbegin(),cigdata.cend(),
		  [&](const pair<char,unsigned> & i)
		  {
		    switch( i.first )
		      {
		      case 'M':
			sum += i.second;
			break;
		      case 'I':
			sum += i.second;
			break;
		      case 'D':
			sum += i.second;
			break;
		      case 'N':
			sum += i.second;
			break;
		      }
		  });
    return sum;
  }

  unsigned insertion_distance( const samrecord & b )
  /*!
    \param b A Sequence::samrecord
    \return The sum of all I elements of a cigar string
  */
  {
    unsigned sum = 0;
    for( samrecord::cigar_iterator i = b.cig_begin() ;
	 i < b.cig_end() ; ++i )
      {
	if ( i->first == 'I')
	  {
	    sum += i->second;
	  }
      }
    return sum;
  }

  unsigned insertion_distance( const bamrecord & b )
  /*!
    \param b A Sequence::bamrecord
    \return The sum of all I elements of a cigar string
  */
  {
    auto cigdata = parse_cigar(b.cigar());
    unsigned sum = 0;
    std::for_each(cigdata.cbegin(),cigdata.cend(),
		  [&](const pair<char,unsigned> & i) {
		    if(i.first == 'I') {
		      sum += i.second;
		    }});
    return sum;
  }

  unsigned deletion_distance( const samrecord & b )
  /*!
    \param b A Sequence::samrecord
    \return The sum of all D elements of a cigar string
  */
  {
    unsigned sum = 0;
    for( samrecord::cigar_iterator i = b.cig_begin() ;
	 i < b.cig_end() ; ++i )
      {
	if ( i->first == 'D')
	  {
	    sum += i->second;
	  }
      }
    return sum;
  }

  unsigned deletion_distance( const bamrecord & b )
  /*!
    \param b A Sequence::bamrecord
    \return The sum of all D elements of a cigar string
  */
  {
    auto cigdata = parse_cigar(b.cigar());
    unsigned sum = 0;
    for_each(cigdata.cbegin(),cigdata.cend(),
	     [&](const pair<char,unsigned> & i) {
	       if ( i.first == 'D') {
		 sum += i.second;
	       }
	     });
    return sum;
  }

  unsigned ngaps( const samrecord & b )
  /*!
    \param b A Sequence::samrecord
    \return Sequence::insertion_distance + Sequence::deletion_distance
  */
  {
    return ( deletion_distance(b) + insertion_distance(b) );
  }

  unsigned ngaps( const bamrecord & b )
  /*!
    \param b A Sequence::bamrecord
    \return Sequence::insertion_distance + Sequence::deletion_distance
  */
  {
    return ( deletion_distance(b) + insertion_distance(b) );
  }

  unsigned mismatches( const samrecord & b )
  /*!
    \param b A Sequence::samrecord
    \return The value of the NM tag of b - Sequence::ngaps.  If the NM 
    tag does not exist, std::numeric_limits<unsigned>::max() is returned.
    Likewise, is Sequence::ngams > the value of the NM tag, 
    std::numeric_limits<unsigned>::max() is returned, as the NM tag
    contains a value that doesn't correspond to the documented definition
    of that tag.

    \note A limitation of the SAM file format is that the number of 
    mismatches are not recorded in any of the required data fields.
    The NM tag is an optional field, and is defined as the edit
    distance of the read from the reference (e.g., the sum of 
    mismatches + indels).  This definition suggests the return
    values described above, but one must trust that the alignment
    software authors have correctly assigned a value to the NM field.
  */
  {
    unsigned sum = 0;
    if ( b.tag_begin() == b.tag_end() ) return numeric_limits<unsigned>::max();

    for( samrecord::tag_iterator i = b.tag_begin() ; i < b.tag_end() ; ++i )
      {
	if ( i->tag() == "NM" )
	  {
	    sum = unsigned (std::stoul( i->value() ) );
	  }
      }
    unsigned ng = ngaps(b);
    if( ng > sum ) return numeric_limits<unsigned>::max();
    return sum - ng;
  }

  unsigned mismatches( const bamrecord & b )
  /*!
    \param b A Sequence::bamrecord
    \return The value of the NM tag of b - Sequence::ngaps.  If the NM 
    tag does not exist, std::numeric_limits<unsigned>::max() is returned.
    Likewise, is Sequence::ngams > the value of the NM tag, 
    std::numeric_limits<unsigned>::max() is returned, as the NM tag
    contains a value that doesn't correspond to the documented definition
    of that tag.

    \note A limitation of the SAM/BAM file format is that the number of 
    mismatches are not recorded in any of the required data fields.
    The NM tag is an optional field, and is defined as the edit
    distance of the read from the reference (e.g., the sum of 
    mismatches + indels).  This definition suggests the return
    values described above, but one must trust that the alignment
    software authors have correctly assigned a value to the NM field.
  */
  {
    bamaux xNM = b.aux("NM");
    unsigned sum = stoul(xNM.value);
    unsigned ng = ngaps(b);
    if( ng > sum ) return numeric_limits<unsigned>::max();
    return sum - ng;
  }
}
