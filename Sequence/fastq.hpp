/*! 
  \file fastq.hpp
  @brief FASTQ class
*/

/*!
  \class Sequence::fastq Sequence/fastq.hpp
  \ingroup seqio
  Publicly derived from Sequence::Seq.
*/
#ifndef __SEQUENCE_FASTQ_HPP__
#define __SEQUENCE_FASTQ_HPP__

#include <Sequence/Seq.hpp>

namespace Sequence
  {
  class fastq : public Seq
    {
    private:
      mutable bool repeat_name;
    public:
      using Seq::Seq;
      std::string quality;
      fastq(void);
      fastq (const std::string &name, const std::string &seq,
	     const std::string & qual);
      fastq (std::string && name, std::string && seq,
	     std::string && qual);
      //! \warning Quality string will be left empty
      fastq (const Seq & s);
      fastq (const fastq & s) = default;
      fastq ( fastq && s) = default;
      //! \warning Quality string will be left empty
      fastq ( Seq && s);
      fastq & operator=(const fastq & ) = default;
      fastq & operator=( fastq && ) = default;
      ~fastq()/*! placeholder for vtable */ {}

      //! Set to true or false for repeating the seq name on third line of output
      void repname(const bool &);
      /*!
	\exception Sequence::SeqException if memory can't be allocated. 
	(This is because the data are temporarily read into char *, 
	because that was found to be faster).
	\exception Sequence::badFormat if the input stream is not
	in FASTQ format
      */
      std::istream & read(std::istream &s);
      /*!
	\param stream a std::ostream
	write the sequence in FASTQ format to \a stream
      */
      std::ostream & print(std::ostream& s) const;
    };
}

#endif

