/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

*/

/*! \file Fasta.hpp
  @brief Declaration of Sequence::Fasta streams
*/

/*!
  \class Sequence::Fasta Sequence/Fasta.hpp
  \ingroup seqio
  Publicly derived from Sequence::Seq, this class defines
  how to read and print sequences in FASTA format, which looks like:\n
  >sequence name 1\n
  ATGATGATCAGATAGACATAGCAGATACATGT\n
  >sequence name 2\n
  ATGTTGGTTTTTTTTTAGAGATGTTTATAGGT\n
  ETC... 
 
  @short FASTA sequence stream
*/

#ifndef FASTA_H
#define FASTA_H

#include <Sequence/Seq.hpp>

namespace Sequence
  {
  class Fasta : public Seq
    {
    private:
    public:
      Fasta():Seq()/*!Generic constructor*/ {}
      Fasta (const std::string &name, const std::string &seq);
      Fasta(const char *name, const char *seq);
      Fasta (const Seq & s);
      ~Fasta()/*! placeholder for vtable */ {}
      /*!
	\exception Sequence::SeqException if memory can't be allocated. 
	(This is because the data are temporarily read into char *, 
	because that was found to be faster).
	\exception Sequence::badFormat if the input stream is not
	in FASTA format
      */
      std::istream&  read(std::istream &s) 
	;
      /*!
	\param stream a std::ostream
	write the sequence in FASTA format to \a stream
      */
      std::ostream& print(std::ostream& s) const;
    };
}
#endif
