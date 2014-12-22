#include <Sequence/Seq8.hpp>
#include <Sequence/util/nibble.hpp>
#include <iostream>

namespace Sequence
{
  Seq8::Seq8(const alphabet_t & _a) : base(),
				      alphabet( _a ),
				      ssize(0)
  {
  }
  
  Seq8::Seq8 ( const std::string & seq,
	       const alphabet_t & _a) : base( std::string(), pack8::dna2vtype(seq,_a) ),
					alphabet( alphabet_t(_a) ),
					ssize(seq.size())
					
  {
  }
  
  Seq8::Seq8( std::string & seq,
	      const alphabet_t & _a) : base( std::string(), pack8::dna2vtype(seq,_a) ),
				       alphabet( alphabet_t(_a) ),
				       ssize(seq.size())
  {
  }
  
  Seq8::Seq8( const std::string & name,
	      const std::string & seq,
	      const alphabet_t & _a) : base( name, pack8::dna2vtype(seq,_a) ),
				       alphabet( alphabet_t(_a) ),
				       ssize(seq.size())
  {
  }
  
  Seq8::Seq8( std::string & name,
	      std::string & seq,
	      const alphabet_t & _a) : base( name, pack8::dna2vtype(seq,_a) ),
				       alphabet( alphabet_t(_a) ),
				       ssize(seq.size())
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
    return ssize;
  }

  std::string::size_type Seq8::length() const
  {
    return ssize;
  }

  std::string Seq8::unpack() const
  {
    return pack8::vtype2dna(std::cref(second),alphabet,ssize);
  }
  
  std::istream & Seq8::read (std::istream & s)
  {
    first.clear();
    second.clear();
    char ch;
    while(!s.eof())
      {
	s.read( &ch, sizeof(char) );
	if( ch != '\0' ) first += ch;
	else break;
      }
    s.read(reinterpret_cast<char*>(&ssize),sizeof( std::string::size_type ));
    second.resize(ssize/2 + (ssize&1));
    s.read(reinterpret_cast<char*>(&second[0]),std::streamsize(second.size()*sizeof(pack8::itype)));
    return s;
  }

  std::ostream & Seq8::print (std::ostream & s) const
  {
    //Write, including the \0
    s.write( first.c_str(), std::streamsize(first.size()+1) );
    s.write( reinterpret_cast<const char *>(&ssize),sizeof( std::string::size_type ));
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
