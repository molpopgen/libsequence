#include <Sequence/util/pack8.hpp>
#include <Sequence/util/nibble.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <functional>
#include <algorithm>
namespace Sequence
{
  namespace pack8
  {
    vtype dna2vtype( const std::string & s, const Sequence::alphabet_t & a )
    {
      using Sequence::nibble::writehi;
      using Sequence::nibble::writelo;
      bool even = (s.size() % 2 == 0.);
      vtype rv(s.size()/2 + !even, 0);
      unsigned j=0;
      for( unsigned i = 0 ; i < s.size() ;  )
	{
	  auto itr = std::find(a.begin(),
			       a.end(),s[i]);
	  if( itr == a.end() )
	    {
	      throw (Sequence::SeqException("Sequence::pack8::dna2vtype error: character not in alphabet"));
	    }
	  alphabet_t::size_type d = 
	    alphabet_t::size_type(std::distance( a.begin(), itr ));
	  writehi(rv[j], itype(d));
	  if( even || i < s.size() - 1 )
	    {
	      itr = std::find(a.begin(),
			      a.end(),s[i+1]);
	      if( itr == a.end() )
		{
		  throw (Sequence::SeqException("Sequence::pack8::dna2vtype error: character not in alphabet"));
		}
	      d = alphabet_t::size_type(std::distance( a.begin(), itr ));
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
    
    vtype dna2vtype( std::string & s, const Sequence::alphabet_t & a )
    {
      vtype rv = dna2vtype(std::cref(s),a);
      s.clear();
      return rv;
    }

    std::string vtype2dna( const vtype & v, const Sequence::alphabet_t & a, const std::string::size_type & len )
    {
      using Sequence::nibble::readhi;
      using Sequence::nibble::readlo;
      std::string rv;
      const bool odd = (len & 1);
      vtype::size_type i = 0;
      const vtype::size_type s =v.size();
      std::for_each(v.begin(),v.end(),
		    [&rv,&a,&i,&s,&odd]( const itype & __i )
		    {
		      rv += a[readhi(__i)];
		      if(i < s-1 )
			{
			  rv += a[readlo(__i)];
			}
		      else if (!odd)
			{
			  rv += a[readlo(__i)];
			}
		      ++i;
		    });
      return rv;
    }
    std::string vtype2dna( vtype & v, const Sequence::alphabet_t & a, const std::string::size_type & len )
    {
      std::string rv = vtype2dna(std::cref(v),a,len);
      v.clear();
      return rv;
    }
  }
}
