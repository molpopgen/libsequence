//! \file Sequence/util/nibble.hpp @brief Declaration of namespace Sequence::nibble */
#ifndef __SEQUENCE_UTIL_NIBBLE_HPP__
#define __SEQUENCE_UTIL_NIBBLE_HPP__

/*! \namespace Sequence::nibble
  @brief Templates for easy writing/reading nibbles
*/
namespace Sequence
{
  namespace nibble 
  {
    /*!
      Write to hi nibble of byte
    */
    template< typename T >
    void writehi( T & byte, T nibble )
    {
      byte = (byte & 0x0F) | T((nibble & 0xF) << 4) ;
    }

    /*!
      Write to lo nibble of byte
    */
    template< typename T >
    void writelo( T & byte, T nibble )
    {
      byte = (byte & 0xF0) | T(nibble & 0xF);
    }

    /*!
      Read high nibble of byte
    */
    template< typename T >
    T readhi( T & byte )
    {
      return (((byte) >> 4) & 0x0F);
    }

    /*!
      Read lo nibble of byte
    */
    template< typename T >
    T readlo( T & byte )
    {
      return (byte & 0x0F);
    }
  }
}


#endif
