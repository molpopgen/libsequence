#ifndef __SEQUENCE_SEQ8_HPP__
#define __SEQUENCE_SEQ8_HPP__

#include <Sequence/Poly8.hpp>
#include <utility>
#include <iosfwd>

namespace Sequence
{
  class Seq8 : public std::pair< std::string, poly8::vtype >
  {
  private:
    std::uint32_t ssize;
  public:
    //! Construct with sequence only                                                                      
    Seq8( const std::string & );
    //! Construct with sequence only.  String will be cleared
    Seq8( std::string & );
    //! Construct with name and sequence
    Seq8( const std::string &,
          const std::string & );
    //! Construct with name and sequence.  Both strings will be cleared
    Seq8( const std::string &,
          const std::string & );
    Seq8(const Seq8 & ) = default;
    Seq8(Seq8 && ) = default;    

    virtual ~Seq8() = default;

    using reference = poly8::vtype::reference;
    using const_reference = poly8::vtype::const_reference;
    using size_type = poly8::vtype::size_type;
    using iterator = poly8::vtype::iterator;
    using const_iterator = poly8::vtype::const_iterator;

    reference & operator[]( const size_type & );
    const_reference & operator[]( const size_type & ) const;
    Seq8 & operator=(const Seq8 & rhs) = default;
    Seq8 & operator=( Seq8 && rhs) = default;
    iterator begin();    
    iterator end();    
    const_iterator begin() const;    
    const_iterator end() const;    
    const_iterator cbegin() const;    
    const_iterator cend() const;    

    std::uint32_t size() const;
    std::uint32_t length() const;

    /*!
      read an object of type Sequence::Seq8 from an istream
    */
    virtual std::istream & read (std::istream & s) = 0;
    /*!
      read an object of type Sequence::Seq8 from an istream
    */
    virtual std::ostream & print (std::ostream & s) const = 0;
  };    

  /*!
    \ingroup operators
    Allows objects derived from Sequence::Seq8
    to be written to output streams.  This operator
    acts by a call to the virtual funtion Sequence::Seq::print
  */
  std::ostream & operator<< (std::ostream & s, const Seq8 & c);
  /*!
    \ingroup operators
    Allows objects derived from Sequence::Seq8
    to be read from output streams.  This operator
    acts by a call to the virtual funtion Sequence::Seq::read
  */
  std::istream & operator>> (std::istream & s, Seq8 &c);
}

#endif
