#include <Sequence/fastq.hpp>
#include <functional>
#include <iostream>

namespace Sequence {
  fastq::fastq(void) : Seq(),quality(std::string()),repeat_name(true)
  {
  }
   
  fastq::fastq (const std::string &name, const std::string &seq,
		const std::string & qual) : Seq(name,seq),quality(qual),repeat_name(true)
  {
  }

  fastq::fastq ( std::string && name, std::string && seq,
		 std::string && qual) : Seq(std::move(name),std::move(seq)),quality(std::move(qual)),repeat_name(true)
  {
  }

  fastq::fastq (const Seq & s) : Seq(s.first,s.second),quality(std::string()),repeat_name(true)
  {
  }

  fastq::fastq ( Seq && s) : Seq(std::move(s)),quality(std::string()),repeat_name(true)
  {
  }

  void fastq::repname(const bool & b)
  {
    repeat_name=b;
  }

  std::istream & fastq::read (std::istream & stream) 
  {
    if( stream.peek() == EOF) return stream;
    if ( char(stream.peek()) != '@' ) 
      throw std::runtime_error("Sequence::fastq::read - error: record did not begin with \'@\'");
    std::string temp;
    stream.ignore(1,'@');
    std::getline(stream,first);
    std::getline(stream,second);
    stream >> std::ws;
    if ( char(stream.peek()) != '+' ) 
      throw std::runtime_error("Sequence::fastq::read - error: third line did not begin with \'+\'");
    stream >> temp >> std::ws;
    if(temp.size() == 1) repeat_name = false;
    quality.resize(second.length());
    stream.read( &quality[0], std::streamsize(second.length()) );
    stream >> std::ws;
    if( second.length() != quality.length() )
      throw std::runtime_error("Sequence::fastq::read - error: sequence and quality strings differ in length");
    return stream;
  }

  std::ostream & fastq::print (std::ostream & stream) const
  {
    stream << '@' << first << '\n'
	   << second << '\n'
	   << '+';
    if( this->repeat_name )
      {
	stream << first;
      }
    stream << '\n'
	   << quality;
    return stream;
  }
} //ns Sequence
