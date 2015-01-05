#include <Sequence/Seq8.hpp>
#include <Sequence/util/nibble.hpp>
#include <iostream>
#include <utility>

namespace Sequence
{
  Seq8::Seq8(const alphabet_t & _a) : base(),
				      alphabet( _a )
  {
  }
  
  Seq8::Seq8 ( const std::string & seq,
	       const alphabet_t & _a) : base( seq.size(), pack8::dna2vtype(seq,_a) ),
					alphabet( alphabet_t(_a) )
					
  {
  }
  
  Seq8::Seq8( unsigned && ssize, pack8::vtype && data, const alphabet_t & _a) : base(std::move(ssize),std::move(data)),
									       alphabet( alphabet_t(_a) )
  {
  }

  Seq8::Seq8( std::string & seq,
	      const alphabet_t & _a) : base( seq.size(), pack8::dna2vtype(seq,_a) ),
				       alphabet( alphabet_t(_a) )
  {
  }

  Seq8::const_reference & Seq8::operator[]( const size_type & i) const
  {
    using nibble::readhi;
    using nibble::readlo;
    return alphabet[ (i&1) ? readlo(second[i/2]) : readhi(second[i/2]) ];
  }

  Seq8::iterator Seq8::begin()
  {
    return second.begin();
  }
  
  Seq8::iterator Seq8::end()
  {
    return second.end();
  }
  
  Seq8::const_iterator Seq8::begin() const
  {
    return second.begin();
  }

  Seq8::const_iterator Seq8::end() const
  {
    return second.end();
  }

  Seq8::const_iterator Seq8::cbegin() const
  {
    return second.cbegin();
  }

  Seq8::const_iterator Seq8::cend() const
  {
    return second.cend();
  }
   
  std::string::size_type Seq8::size() const
  {
    return first;
  }

  std::string::size_type Seq8::length() const
  {
    return first;
  }

  std::string Seq8::unpack() const
  {
    return pack8::vtype2dna(std::cref(second),alphabet,first);
  }
  
  std::istream & Seq8::read (std::istream & s)
  {
    second.clear();
    s.read(reinterpret_cast<char*>(&first),sizeof( std::string::size_type ));
    second.resize(first/2 + (first&1));
    s.read(reinterpret_cast<char*>(&second[0]),std::streamsize(second.size()*sizeof(pack8::itype)));
    return s;
  }

  std::ostream & Seq8::print (std::ostream & s) const
  {
    s.write( reinterpret_cast<const char *>(&first),sizeof( std::string::size_type ));
    s.write( reinterpret_cast<const char *>(&second[0]),std::streamsize(second.size()*sizeof(pack8::itype)) );
    return s;
  }

  std::ostream & operator<< (std::ostream & s, const Seq8 & c)
  {
    return c.print(s);
  }

  std::istream & operator>> (std::istream & s, Seq8 &c)
  {
    return c.read(s);
  }
}
