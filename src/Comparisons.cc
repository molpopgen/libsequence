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
#include <Sequence/SeqExceptions.hpp>
#include <cassert>
#include <string>
#include <cctype>

namespace Sequence
  {
    Mutations TsTv (char i, char j)
    /*!
      Takes two chars, assumed to be nucleotides. The integer returned by this function
      is a member of the enumeration type Sequence::Mutations.
    */
    {
      int k = 0, l = 0;
      switch (char(std::toupper(i)))
        {
        case 'A':
          k = Nucleotides (A);
          break;
        case 'T':
          k = Nucleotides (T);
          break;
        case 'G':
          k = Nucleotides (G);
          break;
        case 'C':
          k = Nucleotides (C);
          break;
        case '-':
          k = Nucleotides (GAP);
          break;
        case 'N':
          k = Nucleotides (N);
          break;
        }
      switch (char(std::toupper(j)))
        {
        case 'A':
          l = Nucleotides (A);
          break;
        case 'T':
          l = Nucleotides (T);
          break;
        case 'G':
          l = Nucleotides (G);
          break;
        case 'C':
          l = Nucleotides (C);
          break;
        case '-':
          l = Nucleotides (GAP);
          break;
        case 'N':
          l = Nucleotides (N);
          break;
        }
      assert(k<=Nucleotides(C) && l <= Nucleotides(C));
      int type = k + l;
      if (type%2 != 0.)	//if odd
        return (Mutations(Tv));	//a transversion
      else if (type%2==0.)	//if even
        return (Mutations(Ts));	//a transition
      return (Mutations(Unknown));	//can be used for error checking
    }

    Mutations TsTv (int i, int j)
    /*!
    Takes two ints, assumed to be integer representations of nucleotides. 
    The way to ensure that the int represents a nucleotide in a valid way is
    to use Sequence::Nucleotides.
    The return value is determined
    by a call to Comparisons::TsTv(int i, int j), where the ints are defined
    in turn by Sequence::Nucleotides
    */
    {
      assert(i<=Nucleotides(C) && j <= Nucleotides(C));
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
		    bool skip_missing,bool nucleic_acid)
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


    unsigned NumDiffs (const std::string & seq1,
		       const std::string & seq2,
		       bool skip_missing
		       ,bool nucleic_acid)
    /*!
      \return the number of differences between two std::strings.  Can skip missing
      data in the same fashion as Comparisons::Different.  If one sequence is shorter
      than the other, the number of positions compared is the length of the shorter 
      sequence.
    */
    {
      unsigned ndiff = 0;
      size_t len = seq1.length();
      if (seq1.length() != seq2.length())
	{
	  len = (seq1.length() < seq2.length()) ? seq1.length() : seq2.length();
	}
      char MISSING = 'N';
      if (!nucleic_acid)
	MISSING = 'X';//for peptide sequences
      for (int i = 0; unsigned (i) < len; ++i)
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
      switch(c)
        {
        case '-':
          return 0;//false, it is a gap...
          break;
        default:
          return 1;
          break;
        }
      return 1;
    }
}

