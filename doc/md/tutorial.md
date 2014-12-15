# Tutorial/overview

[TOC]

This document is a rapid-fire overview of library features.

## Background

I assume a working familiarity with:

* C++ and "C++11"

##Input and output of biological sequences.

A sequence is defined by the "pure virtual" class Sequence::Seq, which publicly inherits from std::pair<std::string,std::string>.  The two members of the pair are the sequence name and the sequence itself, respectively.

A programmer may define new sequences via public inheritance from Sequence::Seq.  The programmer must define the public member funtions Sequence::read and Sequence::print in order to make a valid class.  See Sequence::Fasta and Sequence::fastq for examples.  These two functions (read/print) allow sequences to be read/written to/from C++ streams.

Sequence::Seq defines several functions for data access and various biological operations (Sequence::Seq::Revcom, etc.).

### Reading and writing

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

### Manipulating a sequence

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

#### Gzipped files, etc.

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

## Testing for correct character sets

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

## Input and output of DNA sequence alignments

## Polymorphism tables

## Summary statistics from polymorphism tables

## Coalescent simulation

## Reading "SAM" format records

## Processing BAM files

