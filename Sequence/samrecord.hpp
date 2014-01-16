#ifndef __LIBSEQ_SAMRECORD_HPP__
#define __LIBSEQ_SAMRECORD_HPP__

#include <iosfwd>
#include <string>
#include <vector>
#include <utility>
#include <Sequence/samflag.hpp>

namespace Sequence
{
#ifndef DOXYGEN_SKIP //do not include in library docs
  //fwd declaration of implementation class
  class samrecord_private;
#endif

  struct samtag
  /*!
    \class Sequence::samtag Sequence/samrecord.hpp

    A samtag represents a single TAG:VTYPE:VALUE entry that may appear at the end of a SAM record.

    The struct contains 6 pointers to the beginning and end of the 3 data fields, and member functions
    exist to return those fields as std::string.

    \ingroup HTS
  */
  {
    mutable std::string::const_iterator tag_beg,tag_end,vtype_beg,vtype_end,value_beg,value_end;
    samtag();
    samtag(const std::string::const_iterator & __tag_beg,
	   const std::string::const_iterator & __tag_end,
	   const std::string::const_iterator & __vtype_beg,
	   const std::string::const_iterator & __vtype_end,
	   const std::string::const_iterator & __value_beg,
	   const std::string::const_iterator & __value_end);
    std::string tag() const;
    std::string value() const;
    std::string vtype() const;
    operator std::string();
    operator std::string() const;
  };

  std::ostream & operator<<(std::ostream & o, const samtag & st);

  /*!
    \class Sequence::samrecord Sequence/samrecord.hpp
    \brief A single alignment record from a SAM file

    Intended usage: samtools view (b|s)amfile | ./program_using_this_class

    In essence, this class stores a single line of mapping data as a string.
    This string is then parsed lazily.  Here, lazy parsing means that 
    iterators are stored to the beginning and end of each field, rather than 
    storing each data field as a separate variable.  This drastically speeds
    up the speed of reading in records because, in practice, only a subset of
    the fields may actually be desired.

    When a particular field is wanted, the appropriate member function is called,
    and the iterators are used to construct the desired value.  

    \note The lazy parsing means that repeated calls to retrieve specific fields
    requires constructing a new value from scratch each time.  This is inefficient.  
    So, instead of calling qname() over and over, store the return value in a 
    std::string.

    \ingroup HTS
  */
  class samrecord
  {
  private:
    /*!
      Pointer to implementation class
    */
    samrecord_private * impl;
  public:
    /*!
      Iterator over elements of CIGAR string.  The char is the field type, 
      and the unsigned integer is the value.  Thus, the follwing CIGAR string,
      5M2P5M, may be accessed as follows:

      \code
      Sequence::samrecord s;
      std::cin >> s;
      for( s::cigar_iterator ci = s.cig_begin() ; ci < s.cig_end() ; ++ci )
      {
      std::cout << ci->first << '\t' << ci->second << '\n';
      }
      \endcode
      Will print to the screen:

      M 5

      P 2

      M 5
    */
    typedef std::vector< std::pair<char,unsigned> >::const_iterator cigar_iterator;
    /*!
      Iterator over elements of the optional TAG fields at the end of a SAM record.
      
      Example usage:
      \code
      Sequence::samrecord s;
      std::cin >> s;
      for( s::tag_iterator ti = s.tag_begin() ; ti < s.tag_end() ; ++ti )
      {
      std::cout << ti->tag() << '\t' << ti->vtype() << '\t' << ti->value() << '\n';
      }
      \endcode
    */
    typedef std::vector< samtag >::const_iterator tag_iterator;

    samrecord();
    samrecord(const samrecord & r);
    samrecord(const std::string & s);
    ~samrecord();

    std::string qname() const;
    samflag flag() const;
    std::string rname() const;
    unsigned long pos() const;
    unsigned long mapq() const;
    std::string cigar() const;
    std::string mrnm() const;
    unsigned long mpos() const;
    int isize() const;
    std::string seq() const;
    std::string qual() const;
    std::string tags() const;

    std::istream & read( std::istream & i);
    std::ostream & print( std::ostream & o) const;

    cigar_iterator cig_begin() const;
    cigar_iterator cig_end() const;
    tag_iterator tag_begin() const;
    tag_iterator tag_end() const;
  };

  std::istream &
  operator>>( std::istream & i, samrecord & b );

  std::ostream &
  operator<<( std::ostream & o, const samrecord & b );
}
#endif
