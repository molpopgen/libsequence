#include <Sequence/poly8.hpp>

namespace 
{
  void writehi( Sequence::poly8::itype & byte, Sequence::poly8::itype nibble )
  {
    byte = (byte & 0x0F) | Sequence::poly8::itype((nibble & 0xF) << 4) ;
  }
  
  void writelo( Sequence::poly8::itype & byte, Sequence::poly8::itype nibble )
  {
    byte = (byte & 0xF0) | Sequence::poly8::itype(nibble & 0xF);
  }

  Sequence::poly8::itype emptychar(void)
  {
    Sequence::poly8::itype x = 0;
    writehi(x,Sequence::poly8::itype(Sequence::POLYEOS));
    writelo(x,Sequence::poly8::itype(Sequence::POLYEOS));
    return x;
  }
}

namespace Sequence
{
  namespace poly8
  {
    vtype dna2vtype( const std::string & s )
    {
      vtype rv(s.size()/2 + 1, emptychar());
      bool even = (s.size() % 2 == 0.);
      unsigned j=0;
      for( unsigned i = 0 ; i < s.size() ;  )
	{
	  writehi(rv[j], itype( std::distance( dna_poly_alphabet.begin(),
						std::find(dna_poly_alphabet.begin(),
							  dna_poly_alphabet.end(),s[i]) ) ) );
	  if( even || i < s.size() - 1 )
	    {
	      writelo(rv[j], itype( std::distance( dna_poly_alphabet.begin(),
						   std::find(dna_poly_alphabet.begin(),
							      dna_poly_alphabet.end(),s[i+1]) ) ) );
	      i+=2;
	    }
	  else
	    {
	      ++i;
	    }
	  ++j;
	}
      return rv;
    }

    vtype dna2vtype( std::string & s )
    {
      auto rv = dna2vtype(std::cref(s));
      s.clear();
      return rv;
    }

    vtype dna2vtype( const Sequence::Seq & s ) {
      return dna2vtype(s.second);
    }

    vtype dna2vtype( Sequence::Seq & s ) {
      auto rv = dna2vtype(std::cref(s));
      s.second.clear();
      return rv;
    }

    std::string vtype2dna( const vtype & v ) 
    {
      std::string rv;
      for_each(v.begin(),v.end(),
	       [&rv](const itype & __i) {
		 auto x = dna_poly_alphabet[(((__i) >> 4) & 0x0F)];
		 if(x != '\0')
		   rv += x;
		 x = dna_poly_alphabet[(((__i)) & 0x0F)];
		 if(x != '\0')
		   rv += x;
	       });
      return rv;
    }

    std::string vtype2dna( vtype & v ) 
    {
      std::string rv = vtype2dna(std::cref(v));
      v.clear();
      return rv;
    }
  }
}
