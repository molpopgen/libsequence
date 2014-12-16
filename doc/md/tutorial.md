# Tutorial/overview

[TOC]

This document is a rapid-fire overview of library features.

\section background Background

I assume a working familiarity with:

* C++ and "C++11"

\subsection streams A quick note on streams

In C++, input from streams are handled via operator>> and operator<<. libsequence defines many of these operators for its types (see \ref operators).  These operators are overloaded for uncompressed ASCII streams (aka plain-text files and buffers).  However, all of these operators are compatible with well-designed compressed-file streams, too.  Programmers using libsequence may use [boost](http://www.boost.org)'s filtering_ostream libraries as replacements for <fstream> as they see fit.

\section seq Handling biological sequences.

A sequence is defined by the "pure virtual" class Sequence::Seq, which publicly inherits from std::pair<std::string,std::string>.  The two members of the pair are the sequence name and the sequence itself, respectively.

A programmer may define new sequences via public inheritance from Sequence::Seq.  The programmer must define the public member funtions Sequence::read and Sequence::print in order to make a valid class.  See Sequence::Fasta and Sequence::fastq for examples.  These two functions (read/print) allow sequences to be read/written to/from C++ streams.

Sequence::Seq defines several functions for data access and various biological operations (Sequence::Seq::Revcom, etc.).

\subsection seqio Reading and writing

During read operations, a sequence object may throw an exception of type Sequence::badFormat if the data stream is not in the correct format.  The built-int types (Sequence::Fasta and Sequence::fastq) do not throw on write (although the output stream type could if something bad happens with the stream itself, etc.).

~~~~{.cpp}
#include <Sequence/Fasta.hpp>
#include <fstream>

Sequence::Fasta f;

std::ifstream in("filename.fasta")

try {
in >> f >> std::wd;
} catch (Sequence::badFormat & __b)
{
std::cerr << __b << '\n';
exit(1);
}
~~~~

A sequence may be written to an output stream using the usual C++ output operator <<.

\subsubsection seqgz Gzipped files, etc.

For reading/writing compressed files, see the [boost](http://www.boost.org)'s filtering_ostream libraries.  They have been tested, and "just work" with libsequence objects.

Note: at the time of this writing (and as per my latest testing) the boost libraries do not support opening compressed streams in append mode.  If you need to append to files, then buffer output to a std::ostringstream and write the buffer to a gzFile using [zlib](http://zlib.net) directly:

~~~~{.cpp}
#include <zlib.h>
#include <sstream>
std::ostringstream o;
//Fill o with stuff

gzFile f = gzopen("file.gz","a");
//Write your data to the .gz file
gzwrite(f,o.str().c_str(),o.str().size());
//Clear your buffer
o.str(std::string());
~~~~

\subsection manip Manipulating a sequence

Sequence::Seq provides obvious member functions for element access and iteration. Further, because a sequence inherits from std::pair<std::string,std::string>, a library user may use any of std::string's member functions as well:

~~~~{.cpp}
Sequence fasta f;
//This is the name, whose type is std::string
auto namelen = f.first.size();
//The sequence is also std::string:
auto seqlen = f.second.size();
~~~~

You may also access the data in the sequence (but not the name!) using C++11 range-based for loops:

~~~~~{.cpp}
Sequence::Fasta f("name","TTT");
//The for loop is equavalent to:
//for( auto & c : f.second )
for( auto & c : f )
{
	c = 'C';
}
~~~~~

\section alphabet Testing for correct character sets

The list of characters accepted as DNA is defined in the constant array Sequence::dna_alphabet, declared in Sequence/SeqAlphabets.hpp.  The type of Sequence::dna_alphabet is std::array<const char,17>, and may therefore be iterated over, etc., as per a normal std::array type.

The function Sequence::isDNA accepts a single char as an argument and returns true if the argument is found in Sequence::dna_alphabet:

~~~{.cpp}
#include <Sequence/Fasta.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <algorithm>

Sequence::Fasta f = { "name",	"ATGCZAGC" };  //Z is a non-DNA character

//Find non-DNA characters:
auto itr = std::find_if( f.begin(),f.end(),
     	   		    [](const char & __ch) {
			    return !Sequence::isDNA(__ch);
			   } );

//Delete them from sequences:
  f.second.erase( std::remove_if(f.begin(),
				 f.end(),
				 [](const char & __ch) {
				   return !Sequence::isDNA(__ch);
				 }), f.second.end() );

~~~

The above code block comes from the unit test file alphabets.cc

\section alignments Input and output of DNA sequence alignments

libsequence contains support for I/O and processing of sequence alignments.  The primary use for these functions and classes is to handle the sorts of data that one would collect for population genetic analysis based on Sanger resequencing data and/or phylogenetic analysis.

The relevant modules are:

* Sequence::AlignStream is a template class abstracting the I/O for alignment data.
* namespace Sequence::Alignment

Concrete examples of Sequence::AlignStream include:

* Sequence::ClustalW
* Sequence::phylipData

Developers wishing to handle other input/output formats should study the implementation of those two classes.

As of the time of this writing (2014), these classes may be viewed as quaint.  However, I still think that there is some value to them, and refer the reader to the following sources for usage examples:

* The [analysis](http://github.com/molpopgen/analysis) package that I maintain
* The unit test files AlignStreamTest.cc and AlignmentTest.cc.  These tests provide good coverage of relevant usage cases (and a few cases of user error, too).

\section polymorphism_tables Polymorphism tables

Perhaps the most powerful part of libsequence is the efficient handling of tables of variable sites (polymorphism tables).

There are three basic types for the maninpulation of variation data:

* Sequence::PolyTable is a pure virtual base class for polymorphism tables. In essence, a PolyTable is a std::vector<std::string>, where the strings are the variable sites comprising the haplotypes in the sample.  In addition to this vector of strings, a std::vector<double> stores the positions of the variable sites.  This representation of variation data is the oldest in the library, several functions exist for processing these types.
* Sequence::polymorphicSite is a typedef for std::pair<double, std::string>.  The double represents the site position, and the string represents the state of each individual in the sample.  
* Sequence::Ptable publicly inherits from std::vector<Sequence::polymorphicSite>.  This type is very new to the library, and is the basis of a future rewrite of the code base controlling how variation data are analyzed.

The major difference between Sequence::PolyTable and Sequence::Ptable is how the data are stored internally.  Iteration over a Sequence::PolyTable iterates over _haplotypes_, whereas iteration through a Ptable moves across _variable sites_.  These concepts will become more clear when we look at specific examples below.

\subsection polytable_terms Definitions of terms

Formally, the objects discussed in this section are agnostic with respect to ploidy.  Further, I use the term _haplotype_ here loosely.  If the data that populate a PolyTable or Ptable come from sources such as X-chromosome sequences obtained from males, autosomal sequences from a highly-inbred _Drosophila_ or _Arabidopsis_, or the output of some sort of haplotype phasing algorithm, then the haplotypes are indeed haplotypes (although, for the latter case, one should use the likeliehood of the haplotype inference as a weight on any results, if appropriate).  However, if the input are diploid genotype data, then those data must be split into two strings for that individual (in the case of a PolyTable), which will require arbitrarily assigning the values for a heterozygote to each string.   For such data, __it is user error to then apply any haplotype- or LD-based calculation to the data__.

The only allowed characters in these objects are the set A,G,C,T,N,.,-,0,1.  The first five values should be obvious.  The next two are the identity and gap characters, respectively.  The 0 and 1 may be used in various ways, such as representing arbitrary states of biallelic data, ancestral vs. derived character states, minor/major alleles, or to represent more complex genotypes at a site.  A programmer may check that data contain valid characters using functions declared in Sequence/SeqAlphabets.hpp: Sequence::ambiguousNucleotide and Sequence::invalidPolyChar.

\subsection polytable Sequence::PolyTable in detail

\subsubsection polytables The inheritance hierarchy.

Sequence::PolyTable is a pure virtual class.  As with Sequence::Seq, there are two pure virtuals member functions, Sequence::PolyTable::read and Sequence::PolyTable::print.  A valid class must publicly inherit from Sequence::PolyTable and define these functions.  The library defines the following three classes that publicly inherit from the base class:

* Sequence::SimData is intended to represent binary variation data in the format used by Dick Hudson's coalescent simulation program [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html).  This is the "standard" format used for simulating biallelic sites, and the character states have a very specific meaning: 0 = the ancestral state, 1 = the derived state.  See the example program msstats.cc for how to read these objects in from streams, and the documentation for the file Sequence/SimDataIO.hpp for how to read/write from gzipped streams, binary streams, etc.
* Sequence::PolySites This is a generic/catch-all class for nucleotide-based SNP data.
* Sequence::SimpleSNP This class handles the format used by the software described in [Hudson (2001)](http://www.genetics.org/content/159/4/1805.abstract).

The examples below focus on Sequence::PolySites.

\subsubsection customptables Handling custom input formats

If a user requires reading in data from a custom input format (HapMap, VCF, etc.), take the following steps:

* Create a new class that publicly inherits from Sequence::PolyTable.  For example:

~~~{.cpp}
class myPT : public Sequence::PolyTable
{
};
~~~

* Define the read and print member functions.  See implementations of Sequence::PolySites and/or Sequence::SimData for inspiration.  These functions make your class compatible with operator>> and operator<<, respectively.  Ideally, your class would be compatible with std::istream_iterator<myPT> as well, but that is often a little harder to do.
* Define move-constructors and move-assignment operators.  If you only want to construct myPTs from other myPTs, then the compiler defaults will be fine.  However, if you want to convert from any other type of PolyTable, then you must also define custom move-constructors and move-assigment operators:

~~~{.cpp}
class myPT : public Sequence::PolyTable
{
	myPT( myPT && ) = default;
	myPT( myPT & ) = default;
	myPT( Sequence::PolyTable & ); //you'll need to define this elsewhere
	myPT & operator=( myPT & ) = default;
	myPT & operator=( myPT && ) = default;
	myPT & operator=( Sequence::PolyTable & ); //you'll need to define this elsewhere
	myPT & operator=( Sequence::PolyTable && ); //you'll need to define this elsewhere
};
~~~
* Write unit tests showing that all of the above works. See the following unit tests for examples: testSimDataIO.cc, PolyTableConversions.cc SimpleSNPIO.cc PolySitesIO.cc

An alternative approach to creating your own class is to simply read the data into two vectors:

1. A std::vector<double> representing the mutation positions
2. A std::vector<std::string> representing the haplotypes

From these two containers, a programmer can construct a Sequence::PolySites using move semantics.  The construction of PolyTables is covered in the next section.

\subsubsection polytable_construct Constructing PolyTables

\paragraph polytable_read From streams

Most trivially, a polymorphism table may be read in via Sequence::operator>>, which redirects to Sequence::PolyTable::read.  For example, to read in data from a coalescent simulation that writes to stdout:

~~~{.cpp}
#include <Sequence/SimData.hpp>
#include <iostream>

int main( int argc, char ** argv ) {
    Sequence::SimData d;

    while (! std::cin.eof() )
    {
	std::cin >> d >> std::ws;
	//do something interesting with "d" here:
    }	
}
~~~

See the unit test testSimDataIO.cc for other ways to read/write objects of type Sequence::SimData.

See the unit test PolySitesIO.cc for an example using Sequence::PolySites.

\subsubsection polytable_move Using move construction

Let's say we have the following data in C++ containers:

~~~{.cpp}
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};
  //We can construct a PolySites by moving our data: 
  Sequence::PolySites ps(std::move(pos),std::move(data));

  //Now, accessing the elements of pos and data are undefined!
~~~

These move semantics allow a library user to populate vectors and then fill PolyTables with no extra copying.

Remember that std::move results in the source of the move being left in what may be an undefined state.  Typically, the classes in namespace std will be re-assigned some sensible default.  In the above example, the "pos" and "data" vectors should both be empty.  You can check if that is true on your system by runnning the unit tests in PolyTableTweaking.cc.

In addition to the above example, a well-defined PolyTable has a move constructor and an overload of operator= that takes an lvalue reference as an argument, _e.g._:

~~~{.cpp}
	//The && is an lvalue reference
	myPT & operator=( myPT && ) = default;
~~~

See the header files PolySites.hpp, SimData.hpp, and PolyTable.hpp for examples of declarations.  The corresponding definitions are in PolySites.cc, SimData.cc, and PolyTable.cc, respectively.

\subsubsection polytable_access Accessing the data

The data stored in a PolyTable may be accessed in a variety of ways.  The following types are relevant, and are analogs to the usual typedefs found in STL containers (click to see their documentation):

* Sequence::PolyTable::reference 
* Sequence::PolyTable::const_reference
* Sequence::PolyTable::size_type
* Sequence::PolyTable::data_iterator
* Sequence::PolyTable::const_data_iterator
* Sequence::PolyTable::pos_iterator
* Sequence::PolyTable::const_pos_iterator
* Sequence::PolyTable::const_site_iterator

The following member functions exist for data access:

* Sequence::PolyTable::operator[] returns a reference or const_reference.  The data returned correspond to a std::string representing a haplotype.
* Sequence::PolyTable::begin and Sequence::PolyTable::end return either Sequence::PolyTable::const_data_iterator or Sequence::PolyTabe::data_iterator to haplotypes, depending on the context
* Sequence::PolyTable::cbegin and Sequence::PolyTable::cend return Sequence::PolyTab;e::const_data_iterator to haplotypes
* Sequence::PolyTable::pbegin and Sequence::PolyTable::pend return either pos_iteraor or const_pos_iterator to the mutation positions, depending on the context
* Sequence::PolyTable::pcbegin and Sequence::PolyTable::pcend return const_pos_iterator to the mutation positions.
* Sequence::PolyTable::sbegin and Sequence::PolyTable::send return const_site_iterators
* Sequence::PolyTable::scbegin and Sequence::PolyTable::scend return const_site iterators

The versions with "c" in them may appear redundant, but they are used in C++11 in the context of type deduction using keywords like auto or the declytpe function.

Let's look at some examples.

~~~{.cpp}
Sequence::PolySites p;
//fill p somehow...

//Write the haplotypes to stdout
for( Sequence::PolySites::size_type i = 0 ; i < p.size() ; ++i )
{
	std::cout << p[i] << '\n';
}
~~~

The above example shows us that Sequence::PolyTable::size exists.  There is also Sequence::PolyTable::empty.

Let's do the above using iteration:

~~~{.cpp}
Sequence::PolySites p;
//fill p somehow...

//Write the haplotypes to stdout
for( Sequence::PolySites::data_iterator i = p.begin() ; 
     i < p.end() ; ++i )
{
	//i is an iterator pointing to a std::string
	std::cout << *i << '\n';
}
~~~

Because we have begin and end defined, we can use range-based for loops in C++11:

~~~{.cpp}
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};
  
  Sequence::PolySites ps(std::move(pos),std::move(data)),ps2;

  for( auto d : ps )
  {
	std::cout << d << '\n';
  }
~~~

\subsubsection polytable_manip Manipulating PolyTables

We can take advantage of non-const access to data in order to manipulate what is stored:

~~~{.cpp}
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};

  Sequence::PolySites ps(std::move(pos),std::move(data)),
    ps2(ps);

    //Make everything lower-case
  for( auto & d : ps )
    {
      std::transform(d.begin(),
		     d.end(),
		     d.begin(),
		     [](char & ch) { return std::tolower(ch); });
    }

    //Restore it to upper-case
  for( auto & d : ps )
    {
      std::transform(d.begin(),
		     d.end(),
		     d.begin(),
		     [](char & ch) { return std::toupper(ch); });
    }
~~~

\subsubsection polytable_manip_builtin Methods provided

The library provides a variety of methods for doing things like removing missing data, applying frequency filters, etc.  Unfortunately (for now), these functions are mixed between member functions of Sequence::PolyTable and the file PolyTableFunctions.hpp.  See the documentation for Sequence::PolyTable and PolyTableFunctions.hpp as well as the unit test code PolyTableTweaking.cc for usage examples.  It is possible that a future release of libsequence will deprecate the member functions in favor of standalone functions.

\paragraph polytable_idiot Don't be this guy

Because you have non-const access to the data, you can do this:

~~~{.cpp}
  std::vector<double> pos = {1,2,3,4,5};
  std::vector<std::string> data = {"AAAAA",
				   "AAGAA",
				   "CTGAA",
				   "NAACT"};
  
  Sequence::PolySites ps(std::move(pos),std::move(data)),ps2;

  //You now have a table with unequal haplotype lengths
  //This is bad, and is user error.
  ps[0] = std::string("A");
~~~

\subsubsection polytable_csi Iterating over sites directly

\subsection ptable_detail Sequence::Ptable in detail

\subsection polytable_ptable The relationship between PolyTable and Ptable

\section summstats Summary statistics 

\section coalsim Coalescent simulation

\section sam SAM records

\section bam BAM files

