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
      std::string quality;
      fastq(void);
      fastq (const std::string &name, const std::string &seq,
	     const std::string & qual);
      fastq (const Seq & s);
      fastq (const fastq & s) = default;
      fastq ( fastq && s) = default;
      fastq ( Seq && s);
      fastq & operator=(const fastq & ) = default;
      fastq & operator=( fastq && ) = default;
      ~fastq()/*! placeholder for vtable */ {}
      /*!
	\exception Sequence::SeqException if memory can't be allocated. 
	(This is because the data are temporarily read into char *, 
	because that was found to be faster).
	\exception Sequence::badFormat if the input stream is not
	in FASTQ format
      */
      std::istream & read(std::istream &s) 
	;
      /*!
	\param stream a std::ostream
	write the sequence in FASTQ format to \a stream
      */
      std::ostream & print(std::ostream& s) const;
    };
}

#endif

