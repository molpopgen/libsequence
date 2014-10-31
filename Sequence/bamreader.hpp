#ifdef HAVE_HTSLIB //Will only compile if ./configure detects htslib

#ifndef __SEQUENCE__BAMREADER_HPP__
#define __SEQUENCE__BAMREADER_HPP__

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <Sequence/bamrecord.hpp>

namespace Sequence 
{

  //fwd declaration
  class bamreaderImpl;

  class bamreader
  {
  private:
    std::unique_ptr<bamreaderImpl> __impl;
  public:
    //! Initialize a bamreader object with a file name
    bamreader( const char * bamfilename = nullptr);
    ~bamreader();

    //! \return The next alignment in the file
    bamrecord next_record() const;
    /*!
      \return An alignment record from a specific offset.  Will return an empty record
      if the bgzf_seek functions return an error state.
      \warning If you give an offset that is not the start of a record, the results are undefined
      \note The stream offset is restored to where it was prior to making the call
    */
    bamrecord record_at_pos( std::int64_t ) const;
    //! \return True if bam file is at EOF, false otherwise
    bool eof() const;
    //! \return The return value of bgzf_check_EOF, which checks for whether or not an EOF marker is present in the file.
    bool has_eof() const;
    //! \return True if an error was encountered while reading the bam file, false otherwise
    bool error() const;

    //Stream manipulation

    /*! 
      Rewinds bam file to the beginning (position 0L)
      \return The return value of bgzf_seek
    */
    std::int64_t rewind();
    /*!
      Calls bgzf_seek on the input file
      \return The return value of bgzf_seek
    */
    std::int64_t seek( std::int64_t offset, int whence );
    /*!
      Calls bgzf_close on the input file
      \return The return value of bgzf_close
    */
    int close();

    /*!
      Calls bgzf_tell on input file
      \return The return value of bgzf_tell
    */
    std::int64_t tell();

    //! Iterator type (const only!)
    using refdata_citr = std::vector< std::pair<std::string,std::int32_t> >::const_iterator;
    //! size_type for the container of reference data
    using size_type = std::vector< std::pair<std::string,std::int32_t> >::size_type;
    //! Typedef for reference data (sequence name, length)
    using refdataObj = std::pair<std::string,std::int32_t>;
    //! \return Const iterator pointing to info for first reference sequence
    refdata_citr ref_cbegin() const;
    //! \return Const iterator to end of reference data
    refdata_citr ref_cend() const;
    //! \return the i-th reference name/length pair
    refdataObj operator[](const size_type & i);
    //! \return The complete header
    std::string header() const;
    //! \return The number of sequences in the reference
    std::int32_t n_ref() const;
  };

}

#endif

#endif
