#ifdef HAVE_HTSLIB //Will only compile if ./configure detects htslib

#ifndef __SEQUENCE__BAMRECORD_HPP__
#define __SEQUENCE__BAMRECORD_HPP__

#include <cstdint>
#include <memory>
#include <string>
#include <Sequence/samrecord.hpp>

namespace Sequence
{
  /*! 
    \class Sequence::bamaux Sequence/bamrecord.hpp
    \brief The extra data types at the end of a bam record.
    \ingroup HTS
  */
  struct bamaux 
  {
    /*! 
      The size, in bytes, of value.
      If value_type is a string type,
      the length of the string stored in 
      value will be size + 1 due to \0
      padding at the end.

      \note size will be set to 0
      If the programmer requests the an aux
      data field from and alignment, and that 
      data field does not exist.
   */
    mutable size_t size;
    //! The BAM value_type of the aux data
    mutable char value_type;
    //! The BAM tag type of the aux data
    //mutable std::unique_ptr<char[]> tag;
    mutable char tag[3];
    //! The value of the aux data, in raw bits
    mutable std::string value;
    //Constructors

    /*! For an empty data set.  
      Member data will be set to nonsensical values
      and pointer data set to std::nullptr
    */
    bamaux(); 
    //! For non-empty data
    bamaux( size_t,
	    char[3],
	    char,
	    std::unique_ptr<char[]> & );
    //! Move constructor
    bamaux( bamaux && );
    //! Copy constructor
    //bamaux( const bamaux & );
    //! Assignment operator
    //bamaux & operator=(const bamaux &);
  };

  //!fwd declaration
  class bamrecordImpl;

  /*! 
    \class Sequence::bamrecord Sequence/bamrecord.hpp
    \brief A single alignment record from a BAM file
    \ingroup HTS
  */
  class bamrecord 
  {
  private:
    std::unique_ptr<bamrecordImpl> __impl;
  public:
    //constructors
    bamrecord( std::int32_t blocksize,
	       std::unique_ptr<char[]> && block);
    bamrecord( );
    bamrecord( const bamrecord & );
    //move constructors
    bamrecord(bamrecord&&);
    //destructor
    ~bamrecord();

    bamrecord & operator=(bamrecord&&);
    bamrecord & operator=(const bamrecord&);

    //member functions
    //! True if a record could not be read, or some sort of error occurred
    bool empty() const;

    //Member functions relating to stuff that users care about

    //! \return The read name
    std::string read_name() const;
    //! \return The sequence in base space
    std::string seq() const;
    //! \return Unpacked cigar string
    std::string cigar() const;
    /*!
      \return Quality score string.
      \note Same scale as what is stored in the BAM file
    */
    std::string qual() const;
    //! \return a Sequence::samflag
    samflag flag() const;
    //! \return The mapping quality score
    std::uint32_t mapq() const;
    //! \return The mapping position of the read.  O-offset. -1 = unmapped
    std::int32_t pos() const;
    //! \return The mapping position of the read's mate.  0-offset. -1 = mate unmapped
    std::int32_t next_pos() const;
    //! \return The ID number of the reference sequence where this read maps.  0-offset, -1 = unmapped
    std::int32_t refid() const;
    //! \return The ID number of the reference sequence where this read's mate maps.  0-offset, -1 = mate unmapped
    std::int32_t next_refid() const;
    //! \return Template length
    std::int32_t tlen() const;
    //! \return The length of the read
    std::int32_t l_seq() const;
    //! \return All auxillary data as a single string
    std::string allaux() const;
    //Iterator member functions
    //! Beginning of encoded sequence data. Each integer contains 2 bases, with first base stored in the "high nibble"
    const std::uint8_t * seq_cbegin() const;
    //! One past end of encoded sequence data.
    const std::uint8_t * seq_cend() const;
    //! Beginning of quality data
    const char * qual_cbegin() const;
    //! One past end of quality data
    const char * qual_cend() const;
    //! Returns the record in a raw format
    std::pair< std::int32_t, const char * > raw() const;
    /*! Search auxillary data for a specific tag
      \return The first position of the match if it exists, nullptr if it does not
    */
    const char * hasTag(const char * tag) const;
    /*! Search auxillary data for a specific tag, beginning at position start
      \return The first position of the match if it exists, nullptr if it does not
    */
    bamaux aux(const char * tag) const;
  };

}

#endif

#endif
