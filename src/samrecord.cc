#include <Sequence/samrecord.hpp>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <cstdlib>
/*!
  \defgroup HTS High-throughput sequencing
*/

using namespace std;

namespace Sequence
{
  samtag::samtag()
  {
  }

  samtag::samtag( const std::string::const_iterator & __tag_beg, 
		  const std::string::const_iterator & __tag_end, 
		  const std::string::const_iterator & __vtype_beg,
		  const std::string::const_iterator & __vtype_end,
		  const std::string::const_iterator & __value_beg,
		  const std::string::const_iterator & __value_end) : tag_beg(__tag_beg), 
								     tag_end(__tag_end), 
								     vtype_beg(__vtype_beg),
								     vtype_end(__vtype_end),
								     value_beg(__value_beg),
								     value_end(__value_end)
  {
  }

  string samtag::tag() const
  /*!
    \return name of TAG field
  */
  {
    return string( tag_beg,tag_end );
  }

  string samtag::value() const
  /*!
    \return value of TAG field
  */
  {
    return string( value_beg,value_end );
  }

  string samtag::vtype() const
  /*!
    \return vtype of TAG field
  */
  {
    return string( vtype_beg,vtype_end );
  }

  samtag::operator string()	
  /*
    Typecast to string
  */
  {							    
    return (tag() + ':' + vtype() + ':' + value());		   
  }

  samtag::operator string() const
  /*
    Typecase to string
  */
  {
    return (tag() + ':' + vtype() + ':' + value());
  }

  std::ostream & operator<<(std::ostream & o, const samtag & st)
  /*!
    Writes tag as a string
  */
  {
    o << string(st);
    return o;
  }

#ifndef DOXYGEN_SKIP
  struct samrecord_private
  {
    string record;
    string::const_iterator qname_beg,qname_end,
      rname_beg,rname_end,cigar_beg,cigar_end,
      mrnm_beg,mrnm_end,seq_beg,seq_end,qual_beg,qual_end,
      tags_beg,tags_end,pos_beg,pos_end,mpos_beg,mpos_end,
      mapq_beg,mapq_end,isize_beg,isize_end,flag_beg,flag_end;
    vector< pair<char,unsigned> > cigar_data;
    vector< samtag > samtags;
    char * endptr;
    void parse_record();
    bool valid_cigar();
    void parse_cigar();
    void parse_tags();

    samrecord_private()
    {
    }

    samrecord_private(const std::string & s) : record(s)
    {
      parse_record();
    }
  };

  bool samrecord_private::valid_cigar()
  {
    if ( string(cigar_beg,cigar_end) == "*" )
      {
	return false;
      }
    return true;
  }

  void samrecord_private::parse_record()
  {
    string::const_iterator recend = record.end();
    qname_beg = record.begin();
    qname_end = find_if(qname_beg,recend,::isspace);
    flag_beg = qname_end+1;
    flag_end = find_if(flag_beg+1,recend,::isspace);
    rname_beg = flag_end+1;
    rname_end = find_if(rname_beg+1,recend,::isspace);
    pos_beg = rname_end+1;
    pos_end = find_if(pos_beg+1,recend,::isspace);
    mapq_beg = pos_end+1;
    mapq_end = find_if(mapq_beg+1,recend,::isspace);
    cigar_beg = mapq_end+1;
    cigar_end = find_if(cigar_beg+1,recend,::isspace);
    mrnm_beg = cigar_end+1;
    mrnm_end = find_if(mrnm_beg+1,recend,::isspace);
    mpos_beg = mrnm_end+1;
    mpos_end = find_if(mpos_beg+1,recend,::isspace);
    isize_beg = mpos_end+1;
    isize_end = find_if(isize_beg+1,recend,::isspace);
    seq_beg = isize_end+1;
    seq_end = find_if(seq_beg+1,recend,::isspace);
    qual_beg = seq_end+1;
    qual_end = find_if(qual_beg+1,recend,::isspace);
    tags_beg = (qual_end != recend) ? qual_end+1 : recend;
    tags_end = recend;
    parse_cigar();
    parse_tags();
  }

  void samrecord_private::parse_cigar()
  {
    cigar_data.clear();
    string::const_iterator pibeg = find_if( cigar_beg,cigar_end,::isdigit);
    string::const_iterator piend = find_if( pibeg+1,cigar_end,::isalpha );
    while( pibeg != cigar_end )
      {
	//cigar_data.push_back( make_pair( *piend, atoi( string(pibeg,piend).c_str() ) ) );
	cigar_data.push_back( make_pair( *piend, 
					 strtoul( string(pibeg,piend).c_str(),&endptr,10 ) ) );
	pibeg = find_if( piend+1, cigar_end, ::isdigit );
	piend = find_if( pibeg+1,cigar_end,::isalpha);
      }
  }

  void samrecord_private::parse_tags()
  {
    samtags.clear();
    if(tags_beg==tags_end)return;
    string::const_iterator whitespace_start = tags_beg,
      whitespace_end = find_if( tags_beg,tags_end, ::isspace ),
      c1,c2;

    while( whitespace_start != tags_end )
      {
	c1 = find( whitespace_start, whitespace_end ,':' );
	c2 = find( c1+1, whitespace_end ,':' );
	samtags.push_back( samtag( whitespace_start,c1,
				   c1+1,c2,
				   c2+1,whitespace_end ) );
	whitespace_start = find_if( whitespace_end+1,tags_end, not1(ptr_fun<int,int>(std::isspace)) );
	whitespace_end = find_if( whitespace_start, tags_end, ::isspace );
      }
  }
#endif

  samrecord::samrecord() : impl( new samrecord_private )
  {
  }

  samrecord::samrecord(const samrecord & r) : impl( new samrecord_private(r.impl->record) )
  {
  }

  samrecord::samrecord(const string & s)  
    /*!
      \param s The entire line of a SAM record
    */
    : impl(new samrecord_private(s))
  {
  }

  samrecord::~samrecord()
  {
    delete impl;
  }

  samrecord::cigar_iterator samrecord::cig_begin() const
  /*!
    \return Iterator to beginning of cigar data fields
  */
  {
    return impl->cigar_data.begin();
  }

  samrecord::cigar_iterator samrecord::cig_end() const
  /*!
    \return Iterator to end of cigar data fields
  */
  {
    return impl->cigar_data.end();
  }

  samrecord::tag_iterator samrecord::tag_begin() const
  /*!
    \return Iterator to beginning of parsed tag fields
  */
  {
    return impl->samtags.begin();
  }

  samrecord::tag_iterator samrecord::tag_end() const
  /*!
    \return Iterator to end of parsed tag fields
  */
  {
    return impl->samtags.end();
  }

  std::string samrecord::qname() const
  /*!
    \return QNAME field
  */
  {
    return string(impl->qname_beg,impl->qname_end);
  }

  samflag samrecord::flag() const
  /*!
    \return The SAM flag field as an object of type Sequence::samflag
  */
  {
    return samflag( atoi( string(impl->flag_beg,impl->flag_end).c_str() ) );
  }

  std::string samrecord::rname() const
  /*!
    \return the RNAME field
  */
  {
    return string(impl->rname_beg,impl->rname_end);
  }

  unsigned long samrecord::pos() const
  /*!
    \return the POS field
  */
  {
    return strtoul( string(impl->pos_beg,impl->pos_end).c_str(), &impl->endptr, 10 );
    //return unsigned(atoi(string(impl->pos_beg,impl->pos_end).c_str()));;
  }

  unsigned long samrecord::mapq() const
  /*!
    \return the MAPQ field
  */
  {
    return strtoul(string(impl->mapq_beg,impl->mapq_end).c_str(),&impl->endptr,10);
    //return unsigned(atoi(string(impl->mapq_beg,impl->mapq_end).c_str()));
  }

  std::string samrecord::cigar() const
  /*!
    \return the CIGAR field
  */
  {
    return string(impl->cigar_beg,impl->cigar_end);
  }

  std::string samrecord::mrnm() const
  /*!
    \return the MRNM field
  */
  {
    return string(impl->mrnm_beg,impl->mrnm_end);
  }

  unsigned long samrecord::mpos() const
  /*!
    \return the MPOS field
  */
  {
    return strtoul( string(impl->mpos_beg,impl->mpos_end).c_str(),&impl->endptr,10);
  }

  int samrecord::isize() const
  /*!
    \return the ISIZE field
  */
  {
    return atoi(string(impl->isize_beg,impl->isize_end).c_str());;
  }

  std::string samrecord::seq() const
  /*!
    \return the SEQ field
  */
  {
    return string(impl->seq_beg,impl->seq_end);
  }

  std::string samrecord::qual() const
  /*!
    \return the QUAL field
  */
  {
    return string(impl->qual_beg,impl->qual_end);
  }

  string samrecord::tags() const
  /*!
    \return the TAGS field
  */
  {
    return string(impl->tags_beg,impl->tags_end);
  }

  std::istream & samrecord::read( std::istream & i )
  /*!
    called by operator>>
  */
  {
    getline(i,impl->record);
    impl->parse_record();
    return i;
  }

  std::ostream & samrecord::print( std::ostream & o) const
  /*!
    called by operator<<
  */
  {
    o << impl->record;
    return o;
  }

  std::istream &
  operator>>( std::istream & i, samrecord & b )
  /*!
    \ingroup operators
    Read a Sequence::samrecord
  */
  {
    return b.read(i);
  }

  std::ostream &
  operator<<( std::ostream & o, const samrecord & b )
  /*!
    \ingroup operators
    Write a Sequence::samrecord
  */
  {
    return b.print(o);
  }

}
