#ifndef __SEQUENCE_SEQ8_HPP__
#define __SEQUENCE_SEQ8_HPP__

#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/util/pack8.hpp>
#include <utility>
#include <iosfwd>

namespace Sequence
{
  class Seq8 : public std::pair< unsigned, Sequence::pack8::vtype >
  {
  private:
    alphabet_t alphabet;
    //std::string::size_type ssize;
    using base = std::pair< unsigned, Sequence::pack8::vtype >;
  public:
    //inherit base constructors
    Seq8( const alphabet_t & = dna_alphabet);
    //! Construct with sequence only                                                                      
    Seq8( const std::string &, const alphabet_t & );
    //! Construct with sequence only.  String will be cleared
    Seq8( std::string &, const alphabet_t & );
    Seq8( unsigned && ssize, pack8::vtype && data, const alphabet_t & _a = dna_poly_alphabet);
    Seq8(const Seq8 & ) = default;
    Seq8(Seq8 && ) = default;    

    virtual ~Seq8() = default;

    using reference = alphabet_t::reference;
    using const_reference = alphabet_t::const_reference;
    using size_type = pack8::vtype::size_type;
    using iterator = pack8::vtype::iterator;
    using const_iterator = pack8::vtype::const_iterator;
    using difference_type = pack8::vtype::difference_type;

    /*!
      \return An ASCI character corresponding to the i-th
      position in the unpacked sequence
    */
    const_reference & operator[]( const size_type & i) const;
    Seq8 & operator=(const Seq8 & rhs) = default;
    Seq8 & operator=( Seq8 && rhs) = default;

    /*!
      \return An iterator to the start of the compressed sequence
    */
    iterator begin();
    /*!
      \return The end of the compressed sequence
    */
    iterator end();
    /*!
      \return An iterator to the start of the compressed sequence
    */
    const_iterator begin() const;
    /*!
      \return The end of the compressed sequence
    */
    const_iterator end() const;
    /*!
      \return An iterator to the start of the compressed sequence
    */
    const_iterator cbegin() const;
    /*!
      \return The end of the compressed sequence
    */
    const_iterator cend() const;    

    /*!
      \return The length of the unpacked string, in bytes
    */
    std::string::size_type size() const;
    /*!
      \return The length of the unpacked string, in bytes
    */
    std::string::size_type length() const;

    std::string unpack() const;
    /*!
      read an object of type Sequence::Seq8 from an istream
    */
    virtual std::istream & read (std::istream & s);
    /*!
      read an object of type Sequence::Seq8 from an istream
    */
    virtual std::ostream & print (std::ostream & s) const;
  };    

  /*!
    \ingroup operators
    Allows objects derived from Sequence::Seq8
    to be written to output streams.  This operator
    acts by a call to the virtual funtion Sequence::Seq8::print
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
