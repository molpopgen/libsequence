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

#include <cassert>
#include <algorithm>
#include <Sequence/Seq.hpp>
#include <Sequence/SeqFunctors.hpp>

namespace Sequence
{
  Seq::Seq (void) : SeqBase()
      /*!
	creates an empty sequence
      */
  {
  }
  
  Seq::Seq (const char *name, const char *seq): SeqBase(name,seq)
  {
  }
  
  Seq::Seq (const Seq & seq) : SeqBase(seq.first,seq.second) 
    /*! copy constructor */{}

  std::string Seq::GetName (void) const
    /*!
      Return the sequence name
    */
  {
    return first;
  }

  std::string Seq::GetSeq (void) const
    /*!
      Return the sequence itself
    */
  {
    return second;
  }

  std::string Seq::substr(unsigned beg,unsigned len) const
    /*!
      Mimics the std::string member function of the same name.
    */
  {
    return second.substr(beg,len);
  }

  std::string Seq::substr(unsigned beg) const
    /*!
      Mimics the standardstd::string member function of the same name.
    */
  {
    return second.substr(beg);
  }

  Seq::size_type Seq::length (void) const
    /*!
      Return the total length of the sequence
    */
  {
    return second.length ();
  }

  Seq::reference
  Seq::operator[] (const size_type & i)
    /*!
      Return the i-th element of the sequence.
      \note range-checking is done by assert()
    */
  {
    assert(i < second.length ());
    return second[i];
  }

  Seq::const_reference
  Seq::operator[] (const size_type & i) const
    /*!
      Return the i-th element of the sequence.
      \note range-checking is done by assert()
    */
  {
    assert(i < second.length ());
    return second[i];
  }
  
  bool Seq::operator==(const Seq & rhs) const
  /*!
    \return true if the sequences contain the same data,
    false otherwise.
    \note only the sequences (i.e. this->second and rhs.second) are compared
  */
  {
    return (this->second==rhs.second);
  }

  bool Seq::operator!=(const Seq & rhs) const
  /*!
    \return false if the sequences contain the same data,
    true otherwise.
    \note only the sequences (i.e. this->second and rhs.second) are compared
  */
  {
    return (this->second!=rhs.second);
  }

  Seq::operator std::string() const
  /*!
    allows (implict) cast to std::string
  */
  {
    return second;
  }

  Seq::size_type
  Seq::UngappedLength (void) const
    /*!
      Return length of sequence, excluding the gap character '-'
    */
  {
    size_type ngap = std::count(second.begin(),second.end(),'-');
    return (second.length () - ngap);
  }

  bool Seq::IsGapped (void) const
    /*!
      Returns 1 if the sequence contaings the gap character '-',
      0 otherwise
    */
  {
    return (second.find('-') != std::string::npos);
  }

  void Seq::Subseq (unsigned beg, unsigned length)
    /*!
      \param beg the index along the sequence at which the substring begins
      \param length the length of the subseq
      Acts via std::string.substr().  
      Note that this modifies the data in the object
      by changing thestd::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
      \note range-checking done by assert()
    */
  {
    assert ( beg < second.length() && (beg+length < second.length()));
    second.assign(second.begin()+beg,second.begin()+beg+length);
  }

  void Seq::Complement(void)
    /*!
      Complement the Sequence
      \note This modifies the data in the object
      by changing the std::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
    */
  {
    std::for_each(second.begin(),second.end(),ComplementBase());
  }

  void Seq::Revcom (void)
    /*!
      Reverse and complement the sequence.  
      \note This function modifies the data in the object
      by changing the std::string--if  you want to keep the original sequence, you need to make a copy
      of the object first.
      \return *this
    */
  {
    std::reverse(second.begin(),second.end());
    std::for_each(second.begin(),second.end(),ComplementBase());
  }

  Seq::iterator Seq::begin()
    /*!
      \return an iterator to the beginning of the sequence
    */
  {
    return second.begin();
  }

  Seq::iterator Seq::end()
    /*!
      \return an iterator to the end of the sequence
    */
  {
    return second.end();
  }

  Seq::const_iterator Seq::begin() const
    /*!
      \return a const iterator to the beginning of the sequence
    */
  {
    return second.begin();
  }

  Seq::const_iterator Seq::end() const
    /*!
      \return a const iterator to the end of the sequence
    */
  {
    return second.end();
  }

  const char *Seq::c_str(void) const
    /*!
      \return the the C-style string representing the 
      sequence as a cont char *
    */
  {
    return second.c_str();
  }
}

/*! \mainpage Introduction
  \section purpose Purpose and Intent
  The purpose of this library (which defines namespace Sequence), is to provide
  a set of routines for handling biological sequence data, with an emphasis on how 
  evolutionary geneticists handle data. The intent is not to provide a means of 
  writing sequence-format conversion routines.  In fact, better systems exist for 
  performing those tasks (namely the bioperl project, http://www.bioperl.org).  
  Rather, I intend these libraries 
  to be used as a basis for writing programs for performing many of the computational 
  tasks that are common in evolutionary genetics, a field whose methods are becoming 
  more important to genome analysis and genome comparison.\n
  \n
  Most of the routines are written with nucleotide data in mind, since that is what I
  deal with the most.  The fundamental sequence object is defined by the class 
  Sequence::Seq, which declares a sequence interface and a pure virtual interface for 
  I/O.  There are also routines to translate sequences into peptides
  (Sequence::Translate).\n
  \n
  In practice, sequence data can come in the form of aligned blocks, and the templates
  defined in namespace Sequence::Alignment provides the foundation for dealing with 
  such data.  The virtual base template class Sequence::AlignStream defines an 
  interface for how alignment I/O must work (as Sequence::SeqStream does for
  single sequences).  An example of alignment I/O is defined for ClustalW format alignments
  in the template function Sequence::ClustalW.  
  \n
 
  The library also contains
  contains several classes to do evolutionary genetic analyses.  Classes of
  particular interest are:\n
  1.) Sequence::PolySNP -- analyze molecular population genetic data\n
  2.) Sequence::Comeron95 -- calculate Ka and Ks by Comeron's (1995) scheme\n
  3.) Sequence::Kimura80 -- calculate divergence by Kimura's (1980) method.\n
  \section copyright Copyright and Licensing
  libsequence, copyright Kevin Thornton, University of Chicago, 2002
  \n
  This library is distrubuted under the terms of the GNU public license (GPL) 
  (http://www.gnu.org).
  This means its free, and that you have access to the source code.  And, if you modify the
  library, you must distribute those modifications under the same terms.  The GPL is included
  in the file COPYING in the root of the source directory for the project, and you should read
  it if you have any questions (particularly if you are a commercial user, as it will affect 
  you the most).
  \n
  Most importantly, this library is distributed with no warranty either explicitly stated or 
  implied.
  \section thanks Acknowledgements
  Development of this library had benefited from discussion with several people.  Dick Hudson
  and Eli Stahl provided feedback and much discussion on calculations of summary statistics
  when there are more than 2 states at a site, and Sequence::PolySNP is the result of those
  discussions.  Dick Hudson and Jeff Wall contributed C code that was adapted in to namespace
  Sequence::Recombination.  The coalescent simulation engine is only a slight modification of
  Hudson's original code.
  Gerry Wyckoff provided a table of Grantham's distances that are
  the basis for Sequence::Grantham, and he also provided thousands of comparisons using human/mouse
  divergence to test the output of Sequence::Comeron95.  I should also thank my PhD advisor, Manyuan
  Long, for indulging me the time to work on this when the PCR was running.
  \section portability Compiling the library
  This library has been compiled and tested on a wide variety on Unix systems, including various
  flavors of Linux (http://www.debian.org), Apple's OS X (http://www.apple.com), and Solaris
  systems using g++ 2.9x (http://www.sun.com).  Older versions even compiled under Windows using
  Visual C++, but I don't have access to that platform, and so I will not track portability to it.
  The library is known to compile using gcc 2.9x, 3.x, and 4.x compiler platforms.
 
  As of libsequence 1.5.6, compiler optimizations for Apple G4 and G5 processor systems can be used.
  On a G4, configuring the source code with ./configure --enable-G4=yes sets the options -mcpu=G4,
  -mpowerpc, and -mpowerpc-gpopt.  The option --enable-G5=yes sets -mcpu=G5,-mpowerpc64,-mpowerpc-gpopt

  \subsection compiler Compiler Requirements
  \c libsequence take advantage of many current features of C++, and you compiler needs to support
  them.  Most important amongst these are namespaces, templates (including STL algorithms).
 
  \subsection dependencies Dependencies On Other Libraries
  libsequence requires BOOST (http://www.boost.org) to compile. Note that there are no 
  link-time dependencies on BOOST, only compile-time dependencies.  That means that
  you only need to install the BOOST headers, not the run-time libraries.
 
  \subsection install Installation
  Installing from source is done with the standard 3 commands:\n
  ./configure\n
  make\n
  make install\n
  If you are not familiar with these commands, please consult your local Unix expert.
  \subsection profile Profiling
  Profiling may be enabled by running the configure script with --enable-profile=yes.
  Please remember that accurate profiling of libraries generally requires static linkage (rather than dynamic).
  \subsection debug Debugging
  By default, the library is compiled without debugging symbols, with NDEBUG defined
  (which disables any assertions), and with -O3 to optimize the code.  If you
  wish to enable debugging capabilities, run ./configure with the flag 
  --enable-debug=yes.  The adds -g to the compiler flags, leaves NDEBUG undefined, and 
  does not optimize the resulting object code.  Please note that compiling with debugging
  is only recommended for developers, since it makes the code really big and slow.
 
  \subsection ndebug Notes on NDEBUG
  For those of you unfamiliar with it, NDEBUG is a special symbol for a C/C++ compiler.
  It means "not debugging."  In C, compiling with NDEBUG defined (gcc -DNDEBUG foo.c) disables
  all calls to assert(), and this behavior is identical in C++. By default, the library compiles with -DNDEBUG (see \ref debug).  
 
  \subsection namespaces Namespaces and Scope
  All header files in this library define classes/functions/etc. in namespace Sequence.
  There are also "sub" namespaces, such as Alignment.  None of these are brought into scope
  by default.
  \subsection errors Exceptions and Assertions
  In C++, there are 2 built-in methods to deal with error handling.  The first method is to use
  the assert() function from C, and the second is to use C++ exception handling.  This library
  uses both, but with an emphasis on assertions over exceptions.  The reason for this has to do
  with both efficiency (all the checks to see if we need to throw() an exception can get expensive),
  and code size (including SeqExceptions.h in every file starts to make the library bloated). A 
  better reason, however, has to do with the programming logic.  Much of the code to analyze data
  assumes, for instance, that the data are aligned (implying that all sequences in a data file
  are the same length).  The library provides a function to check if all data read into a vector
  (a vector<Sequence::Fasta *>, for instance) are sequences of the same length (see 
  Sequence::Alignment::IsAlignment).  Thus, it is a programmer error to start analyzing data 
  without first checking that it is aligned, rather than a library error. However, the library
  will check sequence lengths (and a lot of other things), if it is compiled with debugging enabled.
  The checks are done by assert(), and the behavior of assert() is to abort() the program if the
  assertion is false.  Thus, the exceptions thrown by the library deal with errors that a programmer
  cannot reasonably be expected to catch, such as badly formatted data, user input 
  that is unsupported for one reason or another, etc.
 
  \subsection iso ISO C++ Compliance
  As far as I know, everything in this library is up to snuff with respect to ISO C++. All the
  design methods I use are straight from Stoustrup's "The C++ Programming Language" or Meyer's
  "Effective C++" (both from Addison-Wesley).  The library compiles under g++ 3.1.1 
  ( http://gcc.gnu.org ) with both -ansi and -pedantic flags, so that's a good sign at least.
  Reports of any portability problems are appreciated.  Emailing me fixes for the problems
  may actually earn you a beer.  In addition, the coalescent simulation code
  included in this package is implemented in C, and has been modified to successfully compile
  with both -ansi and -pedantic flags.
 
  \subsection thread Thread Safety
  I have never programmed an application using threads, so to be safe, one should assume the
  library is not thread safe.  I will look into this in the future, if I have time.

  \section complink Compiling and Linking Your Code to libsequence
  To compile programs using this library, one must obviously include the appropriate headers
  from the library.  Currently, there is no "lazy man's header" that includes all the headers
  from this package.  The reason for this is discussed in Item 34 of 
  Scott Meyer's book "Effective C++".  Basically, there are a lot of headers, and including them
  all everywhere makes things take forever to compile.
  \n
  To link to the library, use -lsequence -lz when linking up your object code.  The -lz is
  required as of libsequence 1.7.7 as functions in SeqIO.hpp depend on zlib.g/libz.  The
  order of those -l operations does indeed matter!
*/
/*! \defgroup data General I/O
 */
/*!
  \defgroup seqio Classes for sequence I/O
  \ingroup data
*/

/*!
  \namespace Sequence
  The entirety of this library is defined in namespace Sequence.  
  @short The namespace in which this library resides
*/

/*! \page lit References to primary sources
  \section genref General
  Stroustrup, B (1997) The C++ Programming Language, 3rd ed.  Addison-Wesley. -- The standard language reference\n
  Meyers, S (1998) Effective C++, Second Edition. Addison-Wesley -- An essential reference for design using C++
  \section poly Sequence::PolySNP
  Depaulis, F. and M. Veuille (1998) Mol. Biol. Evol. 1788-1790 -- introduction of number of haplotypes and haplotype diversity as tests of neutrality (Sequence::PolySNP::DandVK, Sequence::PolySNP::DandVH)\n
  Fay, JC, and CI Wu (2000) Genetics 155: 1405-1413 -- definition of Sequence::PolySNP::ThetaH\n
  Fu, YX and WH Li (1993) Genetics 693-709 -- definitions of Fu and Li tests\n
  Hudson, R.R., and N. Kaplan (1985) Genetics 111: 147-164 -- definition of Sequence::Recombination::Minrec\n
  Hudson, R.R. (1987) Genetical Research 50:245-250 -- definition of "Hudson's C" (Sequence::PolySNP::HudsonsC)\n
  Hudson, R.R. et al. (1994) Genetics 136:1329-1340 -- discussion of the "Hudson Haplotype test" (Sequence::PolySIM::HudsonsHaplotypeTest)\n
  Simonsen et al.  (1995) Genetics 141: 413 -- power analysis and correction of a Fu and Li test\n
  Tajima, F. (1989) Genetics 123:585-595 -- definition of Tajima's D (Sequence::PolySNP::TajimasD)\n
  Tajima, F. (1993) in Takahata, N. and A.G. Clark (eds) Mechanisms of Molecular Evolution. Sinauer Associates -- general comments on measuring nuceltide diversity (Sequence::PolySNP::ThetaPi), and variances of Pi.\n
  Watterson, G.A. (1975) Theoretical Population Biology 7:256 -- definition of Watterson's Theta (Sequence::PolySNP::ThetaW)\n
 
  \section rec namespace Sequence::Recombination
  Hudson, R.R., and N. Kaplan (1985) Genetics 111: 147-164 -- definition of Sequence::Recombination::Minrec\n
  Hudson, R.R. (1987) Genetical Research 50:245-250 -- definition of "Hudson's C" (Sequence::PolySNP::HudsonsC)\n
  \section com Sequence::Comeron95
  Comeron, J. (1995) J. Molecular Evolution, 41: 1152-1159.\n
 
  \section kim Sequence::Kimura80
  Kimura, M (1980) J. Mol. Evol 16: 111-120.\n
*/
