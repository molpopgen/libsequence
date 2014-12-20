#ifndef __SEQUENCE_UTIL_NIBBLE_HPP__
#define __SEQUENCE_UTIL_NIBBLE_HPP__

namespace Sequence
{
  namespace nibble 
  {
    template< typename T >
    void writehi( T & byte, T nibble )
    {
      byte = (byte & 0x0F) | T((nibble & 0xF) << 4) ;
    }

    template< typename T >
    void writelo( T & byte, T nibble )
    {
      byte = (byte & 0xF0) | T(nibble & 0xF);
    }

    template< typename T >
    T readhi( T & byte )
    {
      return (((byte) >> 4) & 0x0F);
    }

    template< typename T >
    T readlo( T & byte )
    {
      return (byte & 0x0F);
    }
  }
}


#endif
