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

#include <Sequence/SeqEnums.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <string>
#include <cctype>

namespace Sequence
  {
    Mutations TsTv (const char & i,const  char & j)
    /*!
      Takes two chars, assumed to be nucleotides. The integer returned by this function
      is a member of the enumeration type Sequence::Mutations.
    */
    {
      auto k = std::distance( dna_alphabet.begin(),
			      std::find( dna_alphabet.begin(),
					 dna_alphabet.end(),
					 char(std::toupper(i)) ) );
      if ( k > 3 )
	{
	  std::string message("Sequence::TsTv error: ");
	  message += i;
	  message += " is not A,G,C, nor T.";
	  throw SeqException(message.c_str());
	}
      auto l = std::distance( dna_alphabet.begin(),
			      std::find( dna_alphabet.begin(),
					 dna_alphabet.end(),
					 char(std::toupper(j)) ) );
      if ( l > 3 )
	{
	  std::string message("Sequence::TsTv error: ");
	  message += j;
	  message += " is not A,G,C, nor T.";
	  throw SeqException(message.c_str());
	}
      auto type = k + l;
      if (type%2 != 0.)	//if odd
        return (Mutations(Tv));	//a transversion
      else if (type%2==0.)	//if even
        return (Mutations(Ts));	//a transition
      return (Mutations(Unknown));	//can be used for error checking
    }

    Mutations TsTv (const int & i, const int & j)
    /*!
    Takes two ints, assumed to be integer representations of nucleotides. 
    The conversion of int to nucleotide is via Sequence::dna_alphabet
    */
    {
      assert( std::distance(dna_alphabet.begin(),
			    std::find( dna_alphabet.begin(),
				       dna_alphabet.end(),
				       std::toupper(i))) < 4 );
      assert( std::distance(dna_alphabet.begin(),
			    std::find( dna_alphabet.begin(),
				       dna_alphabet.end(),
				       std::toupper(j))) < 4 );
      int type = i + j;
      if (type%2!=0.)	//if odd
        {
          return (Mutations(Tv));	//a transversion
        }
      else if (type%2==0.)	//if even
        {
          return (Mutations(Ts));	//a transition
        }

      return (Mutations(Unknown));	//can be used for error checking
    }

    bool Different (const std::string & seq1,
		    const std::string & seq2,
		    const bool & skip_missing,
		    const bool & nucleic_acid)
    /*!
    Ask if two strings are different.  While this can normally be done by asking
    if (seq1 != seq2) {}, missing data poses a problem here.  If skip-missing == 1,
    missing data (the 'N' character for nucleotide data, 'X' for amino acid)
    are not used to determine if the sequences are different.  If nucleic_acid ==1,
    nucleotide data are assumed, if nucleic_acid==0, protein data are assumed.  
    \note case-insensitive
    \return true if the seqs are different, false otherwise.  If the two sequences
    are of different length, true is returned.
    */
    {
      if(! (seq1.length () == seq2.length ()))
	return true;
      if (skip_missing)
        {
          char MISSING = 'N';
          if (!nucleic_acid)
            MISSING = 'X';//for peptide sequences
          for (unsigned i = 0 ; i < seq1.length() ; ++i)
	    {
	      const char _ch1 = char(std::toupper(seq1[i])),
		_ch2 = char(std::toupper(seq2[i]));
	      if( (_ch1 != MISSING && _ch2 != MISSING) && 
		  _ch1 != _ch2 )
		return 1;
	    }
        }
      else
	{
          for (unsigned i = 0 ; i < seq1.length() ; ++i)
	    {
	      if ( char(std::toupper(seq1[i])) != char(std::toupper(seq2[i])) ) 
		return 1;
	    }
	}
      return 0;
    }


    int NumDiffs (const std::string & seq1,
		  const std::string & seq2,
		  const bool & skip_missing,
		  const bool & nucleic_acid)
    /*!
      \param seq1 A string representing a sequence
      \param seq2 A string representing a sequence
      \param skip_missing If true, missing data characters will not be counted as differences
      \param nucleic_acid.  If true, n/N are the missing data symbol.  If false, x/X.
      \return the number of differences between two std::strings.  Can skip missing
      data in the same fashion as Comparisons::Different.  If one sequence is shorter
      than the other, -1 is returned
    */
    {
      int ndiff = 0;
      size_t len = seq1.length();
      if (seq1.length() != seq2.length())
	{
	  return -1;
	}
      char MISSING = 'N';
      if (!nucleic_acid)
	MISSING = 'X';//for peptide sequences
      for (unsigned i = 0; i < len; ++i)
        {
	  const char _ch1 = char(std::toupper(seq1[i])), 
	    _ch2 = char(std::toupper(seq2[i]));
          if(skip_missing == true)
	    {
	      if( (_ch1 != MISSING && 
		   _ch2 != MISSING) && 
		  (_ch1 != _ch2 ) )
		++ndiff;
	    }
	  else
	    {
              if (_ch1 != _ch2)
                ++ndiff;
	    }
        }
      return ndiff;
    }

    bool Gapped(const std::string &s)
    /*!
    Ask if the std::string contains a gap character.  
    \return true if the string contains gaps, false otherwise
    \note The only gap character checked so far is '-'. Use
    template version for other gap characters
    \deprecated
    */
    {
      return (s.find('-') != std::string::npos);
    }

    bool NotAGap(const char &c)
    /*!
    \return true if a c is not a gap character, false otherwise.
    \note Currently, only '-' is considered to be a gap character
    */
    {
      return c != '-';
    }
}

