# Tutorial/overview

[TOC]

The citation for the library is \cite Thornton:2003vz

This document is a rapid-fire overview of library features.

\section background Background

I assume a working familiarity with:

* C++ and "C++11"

The code snippets below all use C++11 language features (initialization lists, lambda expressions, etc.).  Using GCC/clang, you would need to compile with -std=c++11.  Note, however, that none of the snippets below will actually compile, as they are not complete C++ programs.

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
* Sequence::polySiteVector is a typedef for std::vector<Sequence::polymorphicSite>.  This type has existed for a while, and is a handy way to manipulate data in a SNP-centric way.

The major difference between Sequence::PolyTable and Sequence::polySiteVector is how the data are stored internally.  Iteration over a Sequence::PolyTable iterates over _haplotypes_, whereas iteration through a polySiteVector moves across _variable sites_.  These concepts will become more clear when we look at specific examples below.

\subsection polytable_terms Definitions of terms

Formally, the objects discussed in this section are agnostic with respect to ploidy.  Further, I use the term _haplotype_ here loosely.  If the data that populate a PolyTable or polySiteVector come from sources such as X-chromosome sequences obtained from males, autosomal sequences from a highly-inbred _Drosophila_ or _Arabidopsis_, or the output of some sort of haplotype phasing algorithm, then the haplotypes are indeed haplotypes (although, for the latter case, one should use the likeliehood of the haplotype inference as a weight on any results, if appropriate).  However, if the input are diploid genotype data, then those data must be split into two strings for that individual (in the case of a PolyTable), which will require arbitrarily assigning the values for a heterozygote to each string.   For such data, __it is user error to then apply any haplotype- or LD-based calculation to the data__.

The only allowed characters in these objects are the set A,G,C,T,N,.,-,0,1.  The first five values should be obvious.  The next two are the identity and gap characters, respectively.  The 0 and 1 may be used in various ways, such as representing arbitrary states of biallelic data, ancestral vs. derived character states, minor/major alleles, or to represent more complex genotypes at a site.  A programmer may check that data contain valid characters using functions declared in Sequence/SeqAlphabets.hpp: Sequence::ambiguousNucleotide and Sequence::invalidPolyChar.

\subsection polytable Sequence::PolyTable in detail

\subsubsection polytables The inheritance hierarchy.

Sequence::PolyTable is a pure virtual class that inherits publicly from objects in namespace std:

~~~{.cpp}
//The vector of doubles are the site positions.  The strings are the haplotypes
class Sequence::PolyTable : public std::pair<std::vector<double>, std::vector<std::string> >
{
};
~~~

As with Sequence::Seq, there are two pure virtuals member functions, Sequence::PolyTable::read and Sequence::PolyTable::print.  A valid class must publicly inherit from Sequence::PolyTable and define these functions.  The library defines the following three classes that publicly inherit from the base class:

* Sequence::SimData is intended to represent binary variation data in the format used by Dick Hudson's coalescent simulation program \cite Hudson:2002vy.  This is the "standard" format used for simulating biallelic sites, and the character states have a very specific meaning: 0 = the ancestral state, 1 = the derived state.  See the example program msstats.cc for how to read these objects in from streams, and the documentation for the file Sequence/SimDataIO.hpp for how to read/write from gzipped streams, binary streams, etc.
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
* Sequence::PolyTable::cbegin and Sequence::PolyTable::cend return Sequence::PolyTable::const_data_iterator to haplotypes
* Sequence::PolyTable::pbegin and Sequence::PolyTable::pend return either pos_iteraor or const_pos_iterator to the mutation positions, depending on the context
* Sequence::PolyTable::pcbegin and Sequence::PolyTable::pcend return const_pos_iterator to the mutation positions.
* Sequence::PolyTable::sbegin and Sequence::PolyTable::send return const_site_iterators
* Sequence::PolyTable::scbegin and Sequence::PolyTable::scend return const_site iterators

The versions with "c" in them may appear redundant, but they are used in C++11 in the context of type deduction using keywords like auto or the declytpe function.

__NOTE:__ these iterators are aliases to the underling vectors in the base class, and should be preferred.  In other words, this code:

~~~{.cpp}
Sequence::PolySites ps;

std::for_each(ps.begin(),ps.end(),[](const Sequence::PolySites::const_reference & __s) { std::cout << __s.length() << '\n'; });
~~~

is better than this code:

~~~{.cpp}
Sequence::PolySites ps;

std::for_each(ps.second.begin(),ps.second.end(),[](const Sequence::PolySites::const_reference & __s) { std::cout << __s.length() << '\n'; });
~~~

The reason why will be explained in \ref polytable_csi

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

	//Non-const access to the strings
for( auto & d : ps ) {
	//non-const access to the char within the strings
	for(auto & ch : d ) {
	ch = std::tolower(ch);
	}
}
    //Restore it to upper-case
for( auto & d : ps ) {
	for( auto & ch : d ) {
		ch = std::toupper(ch);
	}
  }
~~~

\subsubsection polytable_manip_builtin Methods provided

The library provides a variety of methods for doing things like removing missing data, applying frequency filters, etc.  Unfortunately (for now), these functions are mixed between member functions of Sequence::PolyTable and the file PolyTableFunctions.hpp.  See the documentation for Sequence::PolyTable and PolyTableFunctions.hpp as well as the unit test code PolyTableTweaking.cc for usage examples.  It is possible that a future release of libsequence will deprecate the member functions in favor of standalone functions.

A user unfamiliar with C++ may think that many features are missing.  How does one permute site positions or the order of the haplotypes?  How can you remove a single haplotype?  These functions are not necessary as they are possible because of the definition of Sequence::PolyTable itself and the functions that already exist in the C++ Standard Template Library (STL).  For example, to remove all haplotypes containing missing data, simply use the erase/remove idiom:

~~~{.cpp}
std::vector<double> pos = {1,2,3,4,5};
std::vector<std::string> data = {"AAAAA",
	"AAGAA",
	"CTGAA",
	"NAACT"}; //This sequence will get erased below

Sequence::PolySites ps(std::move(pos),std::move(data));
ps.second.erase( std::remove_if(ps.begin(),
	ps.end(),
		[](const std::string & __s) {
			return __s.find('N') != std::string::npos;
	}),
	ps.end() );
~~~

Similarly, one may permute the haplotype order and/or the site positions using std::random_shuffle.  If you are not familiar with what is in the STL, it would be a good idea to learn more about it.  I quite like [cppreference](http://en.cppreference.com/w/) as an online source.

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

The examples above iterate over haplotypes.   It is sometimes useful to iterate over the variable sites themselves.  However, iterating over site positions using Sequence::PolyTable::pbegin and Sequence::PolyTable::pend only gives you access to the positions, and not to the character states.  The library allows const access to sites via the member functions Sequence::PolyTable::sbegin and Sequence::PolyTable::send, where the "s" means "site".  These functions return types Sequence::PolyTable::const_site_iterator, whose value_type is __const__ Sequence::polymorphicSite, which is itself a typedef for std::pair< double, std::string >.  Let's look at an example:

~~~{.cpp}
std::vector<double> pos = {1,2,3,4,5};
std::vector<std::string> data =
	{"AAAAA",
	 "AAGAA",
     "CTGAA",
	 "NAACT"};

Sequence::PolySites ps(std::move(pos),std::move(data));

/*
The output of this block will be:
1 AACN
2 AATA
3 AGGA
4 AAAC
5 AAAT
*/
for( auto i = ps.sbegin() ; i < ps.send() ; ++i )
{
std::cout << i->first << ' '
	<< i->second << '\n';
}
~~~

Let's call this type of output "rotated", in the sense that the sites are the rows and the haplotypes are now columns.

The fact that this access is const-only brings up several important points:

1.  Sequence::PolyTable knows if you have accessed it in a const or non-const context.
2.  If you access a PolyTable in a non-const context, the "rotated" representation of the data will be recalculated the next time Sequence::PolyTable::sbegin or Sequence::PolyTable::send is called.  The reason is that your non-const access may have removed some data.
3. In the section called \ref polytable_access, I said that you should prefer iterator-based access to the data via the class's member functions rather than the base classes.  The reason is that _only the former is cabable of distinguishing const from non-const access_, and therefore modifying the data via the latter method will result in const_site_iterators whose data are inconsistent with what the object currently stores.

Let's look at a concrete example of that last point.  The following example appeared above (in the section \ref polytable_manip_builtin):

~~~{.cpp}
std::vector<double> pos = {1,2,3,4,5};
std::vector<std::string> data = {"AAAAA",
	"AAGAA",
	"CTGAA",
	"NAACT"}; //This sequence will get erased below

Sequence::PolySites ps(std::move(pos),std::move(data));
ps.second.erase( std::remove_if(ps.begin(),
	ps.end(),
		[](const std::string & __s) {
			return __s.find('N') != std::string::npos;
	}),
	ps.end() );
/*
The output will be:
1 AAC
2 AAT
3 AGG
4 AAA
5 AAA
*/
std::for_each( ps.sbegin(),ps.send(),
[](const Sequence::polymorphicSite & __p) {
	std::cout << __p.first << ' ' << p.second << '\n';
});
~~~

However, if we had applied the erase/remove step with this code instead:

~~~{.cpp}
/*
Here, we never call Sequence::PolyTable::begin or end.
Rather, we are calling std::vector<std::string>'s versions
of the same functions.
Therefore, does not know that it has been accessed in
a non-const context, and calls to sbegin/send should
be viewed as leading to undefined behavior
*/
ps.second.erase( std::remove_if(ps.second.begin(),
	ps.second.end(),
		[](const std::string & __s) {
			return __s.find('N') != std::string::npos;
	}),
	ps.second.end() );
~~~

How does the above lead to bizarre behavior?  Well, it depends:

1.  If const_site_iterators have never been accessed, then the next call will regenerate the data, and all will be well.
2.  If const_site_iterators _have_ been previously accessed, then the above block will not signal that their data needs to be recalculated.  Exactly what will happen depends on the nature of your non-const access, and should therfore be classified as "scary" at best.

To see how things go badly, look at the unit test file PolyTableBadBehavior.cc

\subsection ptable_detail Sequence::polySiteVector in detail

Sequence::polySiteVector is declared in Sequence/polySiteVector.hpp and defined in polySiteVector.cc.  From these files, you will see that this is an extremely simple class.  Sequence::polySiteVector is simply a std::vector< Sequence::polymorphicSite >, and is therefore related to the objects referred to by Sequence::PolyTable::const_site_iterator (whose value_type is Sequence::polymorphicSite).

This class is most powerful in light of C++11's addition of lambda expressions to the language.  The definition of Sequence::polySiteVector plus the power of lambda expressions leads to a very powerful grammer for manipulating variation tables.  For example, let is remove all sites with invalid characters (as defined in \ref polytable_terms) from a polySiteVector:

~~~{.cpp}
  using psite = Sequence::polymorphicSite;
  Sequence::polySiteVector t = { psite(1.,"AAGC"),
			 psite(2.,"ACZA") }; //site 2 has a non-DNA character

	//This will remove site 2:
  t.erase( std::remove_if( t.begin(),
			t.end(),
			[]( const psite & __p ) {
			     return std::find_if(__p.second.begin(),
						 __p.second.end(),
						 Sequence::invalidPolyChar())
			       != __p.second.end();
			   } ),
	   t.end() );
~~~

The syntax in the above example is compact, readable, efficient, and avoids the pre-C++11 headache of having to define standalone function objects for such simple tasks.  The above code block is from the unit test file polySiteVectorTest.cc.  See the example program polySiteVector_test.cc for some cool usage cases.

\subsection polytable_ptable The relationship between PolyTable and polySiteVector

These two types are intimately-related and may be constructed from one another.

\subsection 

\section summstats Summary statistics

\subsection classic Standard summary statistics
libsequence contains routines for calculating several standard summary statistics (Watterson's \f$\theta\f$, Tajima's \f$\pi\f$, etc.).  The relevant classes are:

* Sequence::PolySNP to calculate statistics from nucleotide data.  The allowed character set is A,G,C,T,N,-.
* Sequence::PolySIM inherits from PolySNP and is intended to be used with biallelic data encoded in a 0/1 format where 0 = ancestral and 1 = the derived character state.  This class is tightly-coupled to Sequence::SimData.  See the example program msstats.cc for how to use this class.

These classes are constructed from objects in the Sequence::PolyTable class hierarchy:

~~~{.cpp}
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolySites.hpp>
#include <iostream>

std::vector<double> pos = {1,2,3,4,5};
std::vector<std::string> data = {"AAAAA",
	"AAGAA",
	"CTGAA",
	"NAACT"}; 

Sequence::PolySites ps(std::move(pos),std::move(data));

Sequence::PolySNP aps(&ps);

//Now, output some summary stats
std::cout << aps.NumPoly() << '\t' //Number of segregating sites
	<< aps.ThetaW() << '\t' //Watterson's theta
	<< aps.ThetaPi() << '\t' //Tajima's pi
	<< aps.TajimasD() << '\n'; //Tajima's D
~~~

Take a look at the class documentation for Sequence::PolySNP and Sequence::PolySIM for a list of all the things that you can calculate from a PolyTable -- there is a lot there.

\subsection classic_fst FST

The class Sequence::FST allows the calculation of \f$F_{st}\f$ statistics from PolyTables + a vector of the sample sizes per population:

~~~{.cpp}
#include <Sequence/PolySites.hpp>
#include <Sequence/FST.hpp>
#include <iostream>
std::vector<double> pos = {1,2,3,4,5};
std::vector<std::string> data = {"AAAAA",
	"AAGAA",
	"CTGAA",
	"CAACT"};
//The sample size is 2 in each subpop:
std::vector<unsigned> sample_sizes = {2,2};
Sequence FST fst_calculator(&data,sample_sizes.size(),&sample_sizes[0]);
std::cout << fst_calculator.HSM() << '\t' //Hudson, Slaktin, Maddison
	<< fst_calculator.HBK() << '\t' //Hudson, Boos, Kaplan
	<< fst_calculator.Slatkin() << '\t'  //Slatkin
	//below are the components of Fst calculations:
	<< fst_calculator.piB() << '\t' //Mean pairwise divergence b/w pops
	<< fst_calculator.piT() << '\t' //Total diversity
	<< fst_calculator.piS() << '\t' //mean within-pop diversity
	<< fst_calculator.piD() << '\n';//The difference between- and within- pop diversity
~~~

\subsection hka The HKA test

The HKA test statistic is available via the functions Sequence::calcHKA.

\subsection stat_future The future 

Sequence::PolySNP and Sequence::PolySIM are "factory" objects, which means that they pre-process your data and contain lots of member functions to calculate various statistics.  A significant problem with this design is that adding new summary statistics breaks the library's compatibility with existing programs compiled against it (because the sizeof(Sequence::PolySNP) changes with the addition of new funtions).   The future of libsequence's interface to summary statistics will be:

* Extracting the preprocessing steps from Sequence::PolySNP to a standalone class
* Rewriting the summary statistic calculations as standalone functions taking the preprocessed data object as a parameter.

The interface described above will be kept because there is a lot of code sitting around that depends upon it.

\section coalsim Coalescent simulation

The sub-namespace Sequence::coalsim contains the routines required for implementing coalescent simulations with recombination using Hudson's algorithm (e.g., the one that underlies his [ms](http://www.ncbi.nlm.nih.gov/pubmed/11847089) program, \cite Hudson:2002vy).  A full introduction to these routines is beyond the scope of this document at the moment, but the namespace has the following features:

* There are no global variables representing the fundamental data structures.  Thus, the code base is prone to fewer side-effects than one would encounter in modifying ms directly.
* It is agnostic with respect to time scale, and may be used for discrete or continuous time scales at the user's discretion
* The current implementation has the Kingman coalescent in mind, in which all coalsecent events are between pairs of lineages.  However, the fundamental data structure (Sequence::coalsim::marginal) will also be compatible with simulating "lamba" coalescents.
* The recombination method implemented in Sequence::coalsim::crossover is Hudson's algorithm.  There is currently no support for the Markovian approximation to this process, but there could be in the future.
* The namespace uses templates to achieve independence from any particular random number generation system.  I have successfully used it with both the C++11 <random> header and the [GSL](http://gnu.org/software/gsl) functions.

The namespace implements several standard/simple demographic scenarios in the file Sequence/Coalescent/DemographicModels.hpp.

The following example programs show more complex use scenarios:

* msmm.cc
* freerec.cc
* fragments.cc
* bottleneck.cc

There is support for simulation involving selection via the header Sequence/Coalescent/Trajectories.hpp.

__DISCLAIMER:__ Please note that this namespace may easily lead to having "too much rope".  As with any simulation interface, knowing how to test what you've coded up is critical, and these functions are intended for people who are comfortable with coalescent theory.

\section hts_tut High-throughput sequencing

\subsection sam SAM records

The class Sequence::samrecord allows processing SAM records from streams:

~~~{.cpp}
#include <Sequence/samrecord.hpp>
#include <iostream>

Sequence::samrecord r;

while( ! std::cin.eof() )
{
	std::cin >> r >> std::ws;
}
~~~

Intended usage for a program using this class would be:

~~~{.sh}
samtools view bamfile | ./program
~~~

The class provides no method for parsing a SAM header.  However, doing so is trivial, and is left to the library user

\subsection samflags SAM flags and bit fields

SAM/BAM data contain "SAM flag" fields.  These fields are 32-bit integers containing a lot of info about the alignment.  They are represented in libsequence by Sequence::samflag.  This type contains boolean variables (Sequence::samflag::is_paired, etc.) representing the various data fields.  The parsing of the bit fields is implemented using data in namespace Sequence::sambits.

\subsection bam BAM files

Sequence::bamreader allows reading directly from BAM files.  The class also supports seeking within a BAM file. An alignment is represented by Sequence::bamrecord, and is returned from a bamreader via Sequence::bamreader::next_record:

~~~{.cpp}
#include <Sequence/bamreader.hpp>

/*
	The header is now parsed if the file was opened successfully
*/
Sequence::bamreader r("file.bam");

if ( r )
{
	while( !r.eof() && !r.error() )
	{
		Sequence::bamrecord rec = r.next_record();
	}
}
~~~

You may access the BAM header info via:

* Sequence::bamreader::header
* Sequence::bamreader::operator[]
* Sequence::bamreader::ref_cbegin()
* Sequence::bamreader::ref_cend()

And use Sequence::bamreader::n_ref to get the number of sequences in the reference.

The Sequence::bamrecord class provides a set of functions to get at the alignment data.  These are direct representations of how the BAM data are stored.  See the class documentation for details.

Some comments:

* Sequence::bamreader is based on the BAM specification.  [htslib](http://htslib.org) is not used for anything other than bgzf decompression and seeking.
* Sequence::bamrecord is move-constructable, meaning that it is lightning fast to copy alignments into containers, etc.

See the author's [pecnv](http://github.com/molpopgen/pecnv) for real-world use of theses classes.  Those programs scan large BAM files in minutes using these classes.
