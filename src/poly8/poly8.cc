#include <Sequence/poly8.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/util/nibble.hpp>
namespace 
{
  Sequence::poly8::itype emptychar(void)
  {
    using Sequence::nibble::writehi;
    using Sequence::nibble::writelo;
    Sequence::poly8::itype x = 0;
    writehi(x,Sequence::poly8::itype(Sequence::POLYEOS));
    writelo(x,Sequence::poly8::itype(Sequence::POLYEOS));
    return x;
  }
}

namespace Sequence
{
  using nibble::writehi;
  using nibble::writelo;
  namespace poly8
  {
    using Sequence::alphabet_t;
    vtype dna2vtype( const std::string & s )
    {
      bool even = (s.size() % 2 == 0.);
      vtype rv(s.size()/2 + !even, emptychar());
      unsigned j=0;
      for( unsigned i = 0 ; i < s.size() ;  )
	{
	  auto itr = std::find(dna_poly_alphabet.begin(),
			       dna_poly_alphabet.end(),s[i]);
	  alphabet_t::size_type d = 
	    alphabet_t::size_type(std::distance( dna_poly_alphabet.begin(), itr ));
	  if( d >= NOTPOLYCHAR )
	    {
	      throw SeqException("character out of range");
	    }
	  writehi(rv[j], itype(d));
	  if( even || i < s.size() - 1 )
	    {
	      itr = std::find(dna_poly_alphabet.begin(),
			      dna_poly_alphabet.end(),s[i+1]);
	      d = alphabet_t::size_type(std::distance( dna_poly_alphabet.begin(), itr ));
	      if( d >= NOTPOLYCHAR )
		{
		  throw SeqException("character out of range");
		}
	      writelo(rv[j], itype(d));
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

    vtype dna2vtype( std::string && s )
    {
      return dna2vtype(std::ref(s));
    }

    vtype dna2vtype( const Sequence::Seq & s ) {
      return dna2vtype(s.second);
    }

    vtype dna2vtype(  Sequence::Seq && s ) {
      return dna2vtype(std::ref(s.second));
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

    std::string vtype2dna( vtype && v ) 
    {
      return vtype2dna(std::ref(v));
    }
  }
}
