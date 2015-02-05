# libsequence revision history

## What this document implies

This document lists major changes between tags/releases of libsequence. The following caveat is important:

If you see a section below referring to libseqeunce version X, and no version X exists as a release, then you must consider the info here inaccurate.  In such cases, the list of changes may refer to planned changes for a future release which may or may not be implemented yet in either the master or some development branch.

## Moving towards a "2.0" release: ABI compatibility in flux

Starting with release 1.8.3, libsequence broke ABI compatibility with previous versions.  The "soname" of the library is now set to 20:0:0, where it will remain __even as future release continue to break the ABI__.   Basically, from here until 2.0, you should consider backwards-compatibility to be in flux.

However:

* As a rule, software depending on libsequence will have to be recompiled, not rewritten.  The ABI breakage is (mostly...) coming from changing the implementation, not the interface.
* If I remove features, they will be marked as deprecated, meaning that they will be removed in 2.0

### TODO:

* Move BAM alphabet to SeqAlphabets
* Change all references of 8 to 4?  Should the name be based on some int8_t or on the fact that it is 4-bit encoded?  8 > 4, so that is obviously a better choice?
* Seq8, PolyTable8 ??
* Explore std::unique_ptr< Sequence::poly8::itype > instead of vector.  It could be faster, but there are implementation/convenience/maintainability costs

1. H12 from Petrov, Messer, et al.
2. G stat for differentiation -- http://arxiv.org/pdf/1403.1552.pdf
3. iHH?
4. http://mbe.oxfordjournals.org/content/31/5/1275.abstract -- DONE
5. Nucleotide, Genotype classes
6. Document: Ptable, bamrecord, bamreader, sam*

### ISSUES:

1. Sequence/SeqRegexes.hpp -- not working.  This will not be fixed until GCC supports <regex>.  The function is now currently implemented in a non-regex manner, which is lame, but it works.

## libsequence 1.8.4

(There is no official 1.8.4 release yet)

* Major improvements to the documentation.  The doxygen output now contains a detailed tutorial.
* fastq stream + fastq unit tests added
* move semantics for the Sequence::Seq hierarchy + unit tests
* Programs in examples/ are now compiled via "make check"
* The build setup for the unit tests is improved
* Old-style enums replaced with C++11 enum classes.  (These enums are WAY better.)
* The enum changes led to many other changes, primarily in the code base dealing wih calculating differences between codons.  The new code is much simpler.  
* A unit test was added to check the calculations of differences bewtween codons.
* Kimura80 and Comeron95 will now return a non-signalling NaN for non-finite values.  Previous versions returned "999", which was a bit silly, in all honesty.
* PolyTable was redefined to publicly inherit from std::pair< std::vector<double>, std::vector<std::string> >.  This change represents a major revision to the code.  This change gives full access to the underlying data, and possibly even "too much rope".  The change also totally wrecks binary compatiblity with previous versions. All unit tests pass.
* Routines for coalescent simulation moved from namespace Sequence to namespace Sequence::coalsim.
* Support for 4-bit encoding of sequence data and variation data
* Sequence::Ptable, introduced in 1.8.3, was a disaster, as it inherited publicly from std::vector, and was therefore a really poor design decision.  The code base has been reverted to a design centered around Sequence::polySiteVector + conversion functions.
* Constructor inheritance for Sequence::Seq enabled.
* Sequence/util/vectorizer.hpp added, in order to facilitate safe inheritance from std::vector
* Issue #6, which was a bug in bamrecord::seq, is fixed

## libsequence 1.8.3 (Dec. 5, 2014)

1.  This release _breaks binary compatibility_ with previous releases.  Programs depending on libsequence will likely need to be recompiled.
2.  New constructor functions were added to Sequence::PolyTable, which breaks binary compatibility.  These new constructors are needed for easier [Rcpp](http://www.rcpp.org/) integration.  Derived classes have been updated to use these new constructors whenever possible.
3.  Some I/O functions are no longer inline.  These include the operator<< and operator>> for Sequence::Seq and Sequence::PolyTable.
4.  The following header files have been removed: Sequence/ensureFloating.hpp (add too much complexity due to metaprogramming), Sequence/preferFloating.hpp (add too much complexity due to metaprogramming), Sequence/RNG/gsl_rng_wrappers.hpp (boost::bind, std::bind, and lambda expressions render these irrelevant.  Removing this header gets rid of the dependency on GSL).
5.  classes Comeron95, PolySNP, and FST no longer inherit from boost::noncopyable.  The classes now declarte copy constructors and the assignment operator= in terms of the new C++11 keyword "delete".  The classes are still not copyable.
6.  All auto_ptr usage has been replaced with unique_ptr 
7.  All function from boost that are now part of C++11 have been replaced with the C++11 function calls.  This mostly includes replacing boost::bind with std::bind and BOOST_STATIC_ASSERT (and the correspponding use of boost/type_traits.hpp) with the C++11 <type_traits> equivalent 
8.  All examples have been updated to C++11, with no dependencies on external libraries like [GSL](http://gnu.org/software/gsl).  Importantly, the GSL is still probably a better way to do random numbers than the C++11 "random" header.  (I personally don't like the design of the C++ random number system).  
9.  Lots of function prototypes have changed. There were some leftover foo(int x) that are now foo(const int & x).  When called many times, the latter may be noticeably faster.
10. A new class, Sequence::Ptable : public std::vector<Sequence::polymorhpicSite> is added.  This is a very powerful class for maninpulating polymorphism data in a general way.
11. The build system has been rebuilt around automake-1.14
12. SimDataIO.hpp/.cc have been rewritten. The new implementation is better, cleaner, and breaks binary compatibility.
13. A new namespace, Sequence::IOhelp, is declared in Sequence/IOhelp.hpp.  It contains some functions to reduce the tedium of some IO conversions.
14. A new, faster implementation of Sequence::SimData::read
15. The configure script how checks for the existence of [htslib](http://www.htslib.org)
16. Sequence/bamreader.hpp and Sequence/bamrecord.hpp added for direct reading from bam files.  Will only compile if [htslib](http://htslib.org) is detected by the configure script
17. I've performed extensive code auditing using cppcheck.  Many issues were resolved.  Most were stylistic, but they do make the code more idiomatic.
18. A unit test suite is starting to take shape in the test/ subdirectory.  These test programs depend on the [boost](http://www.boost.org) "test" library.  Future tests could add future dependencies.
19. The unit testing revealed a variety of subtle issues in the code base. All of these were in "dark corners" of the library, _e.g._ functions that have seen little or no use in production programs.  

## Historical data (copied from the readme.txt previously found at molpopgen.org)

current version 1.8.0 (Jan 14, 2014).  Two files needed #include <sstream> to used std::ostringstream.  Previous library versions were failing to compile on some systems due to the missing header.

version 1.7.9 (Dec 12, 2013)  Forward declaration of templates in namespace std were removed from Sequence/typedefs.hpp.  The library will now compile on OS X Mavericks.

version 1.7.8 (Nov 8, 2013):  Fixed (a possibly 10 year old)  bug in SingleSub.cc.  This bug misclassified AGG to CGG changes as q_0 rather than q_2V.  This would lead to errors in synonymous/nonsynonymous substitution counting.  Programs like gestimator and polydNdS were affected.  However, I have run through a lot of my own data, and recorded no differences in output.  You have to have had this specific change (or possibly others involving amino acids with 6 codons in the family).  The SingleSub code is now vastly simplified and gives correct output (fingers crossed).

Binary compatilibity with previous versions of libsequence is now broken:   All exception specifications are removed in anticipation of that language feature being removed/deprecated. 

A new file, SimDataIO.hpp has been added, allowing .gz and binary I/O for SimData objects.  See examples/test_SimDataIO.cc for documentation.  This now changes how one compiles programs depending on libsequence:  you need zlib installed and must include -lsequence -lz (in that order) when compiling programs.

version 1.7.7: never happened. never publicly released

version 1.7.6 (Feb 13, 2013)  Fixed inconsistent exception declarations in templates.  libsequence is now compilable with the clang compiler.  To use: CXX=clang++ CC=clang ./configure; make ; sudo make install

version 1.7.5 (Dec. 18 2012).  The operator>> for Clustalw was not processing whitespace properly when compiled with new versions of gcc.  That issue is now fixed.

version 1.7.4 (Oct 31, 2011).  Replaced all instances of typename ptrdiff_t with std::ptrdiff_t, as gcc 4.6 is picky enough to care.

1.7.3 (Aug 11, 2011).  The more efficient version of Sequence::Recombination::Disequilibrium failed to set a flag in the return value, leading to duplicated output.  Now fixed.

version 1.7.2 (Aug 11, 2011). The bugfix in 1.7.1 didn't really work.  Implementation class data weren't fully-copied.  That is now fixed

version 1.7.1 (Aug 3, 2011).  Fixed bug in samrecord class.  Copying of samrecords led to segfaults in previous versions, due to the private pointer-to-implementation-class not getting its data copied to the new object.

version 1.7.0 (December 13, 2010).  A bug is fixed in the calculation of the number of states at a SNP when doing LD calculations.  For monomorphic sites, the old version had problems.  Also, an initial release of sam/bam record parsing is included (Sequence/sam*.hpp)

version 1.6.9 (June 25,2010) A bug is fixed in calculation of haplotype diversity.  When there is considerable missing data in a sequence file, this statistic was severely under-estimated.  Previous analyses using this statistic on real data should be re run.

version 1.6.8 (March 29, 2009)  Implemented a much more efficient means of calculating basic LD statistics.  This is not (yet) available through the PolySNP/PolySIM classes yet, but must be accessed via the header file Sequence/Recombination.hpp

version 1.6.7 (March 24, 2009)  Implemented a safer iteration through coalescent histories, in light of a recent bug in "ms", on which my routine was based.  Also implemented a "max marker distance" for LD calculations.  NOTE:  This version is binary-incompatible with previous releases, meaning any software depending on libsequence should be recompiled.

version 1.6.6 (Feb 4, 2009)  Updated to be compile cleanly with gcc 4.3.x.  

version 1.6.5 (Nov. 29, 2007) Updated math function calls to C99 standard.  Should compile find on Apple's "Leopard" operating system now.

version 1.6.4 (Jun 19, 2007). Fixed bug where haplotype statistics for SNP data were incorrectly calculated in some cases. When alignment columns containing polymorphic sites also contain a mix of upper- and lower-case letters, the extra case was incorrectly detected as an extra haplotype. This affected the following statistics: # haplotypes, haplotype diversity, Wall's statistics, and Hudson's (1987) estimator of 4Nr.

version 1.6.3 (Jun 22, 2006).  Several changes.  Templates for
coalescent simulation should now be much more flexible with respect to
random number generators.  Also, the new design allows them to work with
boost::bind + GSL routines, simplifying life quite a bit.  Code for Hudson's
Snn test has been added, as well as code for simulating Coop &amp;
Griffiths-style trajectories of beneficial alleles.

version 1.6.2 (March 5, 2006).  Three bugs fixed.  First, in the FST module, when there are lots of populations, the sum of all the population weights was sometimes detected as not summing to 1, when they actually did.  The bug was not accounting for numerical precision, which is now fixed.  Second, a rare segfault in the sliding window module (PolyTableSlice.hpp) is now fixed.  Third, pairwise LD statistics were being calculated incorrectly in the presence of an outgroup.  Now fixed.

version 1.6.1 (Jan 6, 2006).  New summary statistics added for SNP data.  
Some example demographic models for coalescent simulations added (see header 
Sequence/Coalescent/DemographicModels.hpp).  Five bugs fixed:
  1.) src/CoalescentInitialize.cc: Fixed bug in init_sample. Regions with 0 "sites" now have the last "site" in the chromosome labelled 0, which is correct. Previously, it was labelled 1. 

  2.) src/HKA.cc: Fixed minor error in calculation of f_hat which was pointed out by Andy Kern. I was taking the mean(f_hat) across loci as the estimate, rather than following equation 5 of the HKA paper closely enough. The effect was small, and estiamtes of theta and T_hat were unaffected or affected only very little, repsectively.

  3.) src/PolySNP.cc: Fixed glitch in calculation of ThetaW. When totMuts==false, only the number of bi-allelic SNPs was used in the numerator of ThetaW. Now, the total number of variable sites is used, as intended.

  4.) src/Recombination.cc (affects calculation of LD statistics): In the calculation of D', the 11 gamete was ancestral/ancestral, rather than derived/derived. This is now fixed to match the documentation.

  5.) src/CoalescentTreeOperations.cc: in function total_time_on_arg, there was a bug leading to an exception always being thrown. This was due to an inequality being tested as &lt; , when &gt; was what was intended.

version 1.6.0 (Oct 15, 2005).  In class PolySNP, summary statistics
that require outgroup data now return "nan" when the outgroup is not
present.  Bug fix: when a PolySNP/PolySIM object was constructed from a
SimData object that has 0 segregating sites, the value of 0 was returned for
number of haplotypes.  It now correctly returns 1.

version 1.5.9 (Sept 21, 2005)  more efficient placing of mutations
on simulated genealogies.  class SimData can now be used as a base class.
Note: binary compatibility broken w/previous versions!

version 1.5.8 (July 21, 2005) fixed memory leak introduced in 1.5.7

version 1.5.7 (July 19, 2005)  more efficient crossover routine for
coalescent simulation.  Bug fixed in Sequence::SimpleSNP that lead to some
columns mistakenly getting dropped on input.

version 1.5.6 (Jun 15, 2005).  No bugfixes.  A function in namespace
Sequence::Alignment is no longer recursive.  Library source can now be
configured to generate cpu-specific optimizations on G4 and G5 powerpc
systems running OS X with gcc.

version 1.5.5 (Jun 9, 2005).  Bugfix in Sequence::PolySNP where
preprocessing was not correctly applied.  The effect of the bug was that
some results would be incorrect when performed in certain orders.  This
affected no programs that I have written and distribute, and analysis of
simulated data (using Sequence::PolySIM) was unaffected.  In short, it
shouldn't have hurt anybody.  The non-const operator[] of
Sequence::PolyTable now sets non_const_access to true (this fixes a minor
bug, which again will rarely affect anything because most operations in the
library are on const PolyTables).  



version 1.5.4 (May 27, 2005).  Major new additions.  Old coalescent
simulation code is removed.  New types and functions allowing simulations
with recombination added.  changes in support for random number generators
implemented (boost-based stu
ff removed in favor of the GNU scientific
library).  This release breaks runtime compatibility.


version 1.5.3 (April 21, 2005) bug in PolySNP::ThetaPi() (i.e.
calculation of pi for sequence data) where sites in tables where all but 1
individual had missing data were not handled correctly

version 1.5.2 (April 8, 2005) The bugfix in 1.5.1 was incomplete.
Now fixed and tested. 

version 1.5.1 (April 5, 2005).  Bugfix in sliding window code for
non-overlapping windows along the DNA sequence.  New example code.  Cleanup
to code base so that .cc files compile without preprocessor macros.

version 1.5.0 (Jan 11, 2005). Several minor bugfixes: implicit typenames fixed so that library now compiles under gcc 3.4, a bug in the documentation of the formula for PolySNP::ThetaPi() is now corrected, &lt;ieeefp.h&gt; is now included in the HKA source code for systems that have that header (i.e. cygwin), and the estimate of the sample variance in &lt;Sequence/descriptiveStats.hpp&gt; is now unbiased (i.e uses n-1 rather than n in the denominator). In addition, there is a policy change in the calculation of several summary statistics in the classes Sequence::PolySNP and Sequence::PolySIM. When a summary statistic is undefined, the value nan (not a number) is now returned.

version 1.4.9 (Dec 3, 2004).  Minor bug in HKA code node fixed.
That probably affected nobody, as I don't distribute any programs that
depend on it.  In addition, several issues with compiling under gcc 3.4 were
fixed.  These were not related to algorithms, etc., but rather involed
non-standard declarations and such.

version 1.4.8 (Nov 3, 2004) -- Major bug in counting derived states
fixed in Sequence::PolySNP.  Analyses using previous versions of the library
should be repeated.  Analysis of coalescent sim data (using
Sequence::PolySIM) was not affected.


version 1.4.7 (Oct 7, 2004) Class Sequence::PolySIM now correctly
sets both Wall's B and Q to -1 if there are &lt; 2 segregating sites.
Previously, only Wall's B was set correctly

version 1.4.6 (Oct 6, 2004)
Fixed bug in HKA calculations where the poly data in species 1 had 0 seg
sites.  Sequence::SimpleSNP now does better checking of its input.
PolyTableSlice (aka "sliding windows") now stores windows both with and
without segregating sites.

version 1.4.5 (Aug 22, 2004).
Bugs fixed:
The output shift operator for template class AlignStream now takes a const reference.    
The input routines for SimParams now reserve memory, which seems to have  
fixed a rare crash on OS X systems.  

New features:
Sequence/HKA.hpp provides an interface to calculating the HKA test of
Hudson, Kreitman, and Aguade.

The class Sequence::PolySIM (Sequence/PolySIM.hpp) is now much more
efficient, due to a reorganization of code in the base class
Sequence::PolySNP (Sequence/PolySNP.hpp).  In some applications,
calculations can be up to 1/3 faster.  The reorganization of code breaks
binary compatibility with previous library versions.

Sequence::SimParams (Sequence/SimParams.hpp) now has a member function
fromfile(FILE *), mimicing that of Sequence::SimData.

version 1.4.4 (Jul 8, 2004) : minor bugfix.  All files that use
toupper now include &lt;cctype&gt;

version 1.4.3 (Jun 4, 2004) :
   * bugfixes : in PolySNP.cc, NumSingletons and NumExternalMutations could
give incorrect results with missing data (in cases were all but 1 individual
is missing).  This affected the output of Fu &amp; Li statistics.  
   * class SimData now does output exactly like hudson's ms program when
there are 0 segregating sites
   * class SimData can now be read in from C-style FILE *.
      * made implementations of routines for coalescent simulation a little
more robust to uniform random number generators that return 1.0
   * in Alignment.hpp, an exception specification in a declaration did not
match the definition.  Fixed.

version 1.4.2 (May 14, 2004) -- several bug fixes &amp; new features:
   *removal of unused variables in several place
      *efficiency improvement in Sequence::NumDiffs() (in Comparisons.hpp)
         *updated Sequence::uniform (in BoostRandomNumbers.hpp) so that it
compiles and is compatible with std::random_shuffle
   *new functions added to namespace Alignment and PolyTableFunctions.hpp
      *added several new headers providing template functions for simple
statistics
   *fixed some bugs in Sequence::SimpleSNP (sample size now printed
correctly upon output when there is an outgroup)
   *fixed bug in Recombination::Disequilibrium (frequency filter now works,
case-insensitivity implemented).  And, for D and D', the 11 gamete type is
in terms of the minor (or derived) allele
   *fixed use of assert() to deal with compiler warnings

version 1.4.1 (Apr 26, 2004) -- Many bug fixes!
   * PolySNP -- member functions explicitly check for gaps in snps
      * PolySNP -- NumExternalMutations() was incorrect, leading to bad
values for Fu/Li statistics with outgroup.  Fixed.
   * PolySNP/PolySIM -- bugs fixed in calculations of Wall's statistics.  Q
was simply incorrect, and all statistics were incorrect for polymorphism
tables constructed with invariant ingroup data
   * SimpleSNP -- fixed segfault on input for files containing outgroup data
      * PolyTableFunctions.hpp added, contains functions to manipulate
PolyTable objects in non-const ways

version 1.4.0 (Apr 16, 2004) -- 
   * PolySNP::FuLiD() is now correct, no longer returning -inf
      * PolySites::read() is now public
         * PolyFunctional.hpp is now included
	 

version 1.3.9 (Apr 12, 2004) -- version 1.3.8 had a bug where the
Translate function printed all the codons it processed to stdout.  That's
fixed in this version.  Also fixed a design flaw where output for PolyTables
did not work in const contexts.

version 1.3.8 (Apr 08, 2004) : 2 bug fixes: The first was in
Sequence::PolyTable::Binary(), and resulted in the binary-format data not
being assigned to the object. The second bug was in Sequence::SimData, and
resulted in "ms output" being read in incorrectly when using C++-style
input. The member function SimData::fromstdin(), which uses C-style input
for speed, however, worked just fine.  Other new features include : 
   * non-const operator[] added to PolyTable
   * public memeber assign() added to AlignStream
   * print() defined for Clustalw
   * removed Hudson's code (since in practice everyone just pipes
the results of ms to another program anyways...)
   * libsequence should now be case-preserving of input data
   * PolyTable now has assign() member functions, to mimic std::vector
   * PolyTableSlice, a template class for sliding window analyses, is
completely re-written and easier to use
   * updates to documentation
   

version 1.3.7 (Feb 11, 2004) : Sequence::SimData has 2 changes.  First: you
can now say while(cin &gt;&gt; d) {} and it works as expected, where d is an
instance of SimData.  Two, SimData::fromstdin() now returns the return value
from the C function fscanf, which makes checking for end of input more robust.
Also, the class Sequence::uniform in &lt;Sequence/BoostRandomNumbers.hpp&gt; had
it's declaration changed to be compatible with old and new verions of boost.

version 1.3.6 (Feb 02, 2004) : new constructor for PolyTable,
allowing construction from a range of PolyTable::const_site_iterator.  
New classes in &lt;Sequence/shortestPath.hpp&gt; to aid in dealing w/differences
between codons.

version 1.3.5 (Jan 19, 2004) -- Bugfix: PolySNP::Minrec() was
returning a value for the Rmin statistic that was incorrect for a small
number of data sets.  This has been fixed.  Other changes : random number
system added, based on boost classes.  Autoconf macros now work properly at
creating #define's for needed headers.

version 1.3.4 (Jan 12, 2004) -- Binary compatibility broken w/previous version.
PolyTable column iterators added.  New function PolyTable::RemoveAmbiguous
added. Ka/Ks code cleaned up.  See ChangeLog and manual for details.
Sequence::Seq now inherits from std::pair&lt;std::string,std::string&gt;, and
various template functions were updated to accomodate that.  Template
memebers added to Sequence::PolyTable.
Requires boost (http://www.boost.org) to compile.

version 1.3.3 (Nov 1, 2003) -- Exception specifications added.
Development headers now included in subdirectory Sequence (so you now need
to say things like (#include &lt;Sequence/Fasta.hpp&gt;).

version 1.3.2 (August 7, 2003) -- Yesterday's release forgot to
include the iterator types for Sequence::AlignStream&lt;T&gt;.  Here they are.  
The soname of the shared library is unchanged, since no one downloaded 1.3.1
yet :).

version 1.3.1 (August 6, 2003) -- updated to autoconf 2.5.  added
new iterator type to Sequence::PolyTable.  soname changed to 12:0:0.
Sequence::FST updated.

version 1.3.0 (July 9, 2003) -- Bugfix for
Sequence::Comparisons::NumDiffs().  Function now skips missing data as it
should. The bug was due to a rewrite of the fxn several versions ago.

version 1.2.9 (Jun 21, 2003) -- Bugfix for PolySNP::ThetaW().  Caused
the denominator to be incorrectly calculated.  It only occured when running
my 'compute' program (from analysis 0.4.8) on several files at once. May
have been gcc-3.3 specific, but now it works more generally.

version 1.2.8 (Jun 1, 2003) -- no major changes, so no particular
need to update to this version.  Restores compatibility with gcc 2.95

version 1.2.7 (Jun 1, 2003) -- added label(unsigned) to
Sequence::Hudson2001.  Also fixed a preprocessor bug that prevented
SeqUtilitites.hpp from being included. No programs were affected by that
bug.  Also cleaned up example code.

version 1.2.6 (May 21, 2003) -- Two things.  One, a bug in PolySNP::ThetaH is fixed which
resulted in the value of the statistic being about 6 times too large.
(Calculation of ThetaH for simulated was fine). Second, several substantial
speed improvements were made to class PolySNP, which speed up the analysis
of sequence data.  First, calculation of Rmin is no longer recursive, which
led to a 10-fold speed up on my machine.  Second, a preprocessing of the
frequency spectrum now allows summary stastics to be calculated in linear
time (at the cost of a bit more memory).

version 1.2.5 -- fixed missing #include&lt;cassert&gt; in Comeron95.cc

version 1.2.4 -- Several bugfixes.  PolySNP::WallStats went beyond
array bounds for the number of segsites &lt; 2.  This is now fixed, and -1 is
returned for the statistics, indicating not applicable.  Also, several
missing #include directives were added to source files, which should make
compilers more happy (don't know why gcc 3.3 wasn't complaining, but 2.95
did complain).

version 1.2.3 -- Bugfix release for OS X systems.
Sequence::PolyTable::Assign(double *, unsigned, char **, unsigned) ran into
an issue where the c++ compiler on OS X does not seem to generate
standard-compliant code for vector&lt;T&gt;::assign().  A workaround was
implemented.


version 1.2.2 -- Template specializations of routines in namespace
Sequence::Alignment are now provided for std::string.  A bug causing a
segfault was fixed in Sequence::Alignment::RemoveGaps&lt;T&gt;()

version 1.2.1 -- loads of changes. Sequence::Poly is removed,
replaced by Sequence::PolySNP and Sequence::PolySIM.  Sequence::FST is
added.  see ChangeLog for details

version 1.2.0 -- lots of changes, mostly internal to the code.  See ChangeLog.

version 1.1.9 -- lots of changes.  Iterators added to Seq and
PolyTable.  See ChangeLog for details.

version 1.1.8 -- Minimum # of recombination event and disequilibrium
coefficients now work with sequence data.

version 1.1.7 -- many changes, see ChangeLog for Details

version 1.1.6 -- fixed bug in input routine for Hudson2001,
resulting in occasional segfaults

version 1.1.5 -- fixed bugs in
PolyTable::RemoveMultiHits(),ApplyFreqFilter(), and RemoveMissing() that led
to segfaults.  Also, members of Poly that returned int now return unsigned int

version 1.1.4 -- general code cleanup.  Also, removed a
simplification in the calculation of Sequence::Poly::ThetaW(SimData *)

version 1.1.3 -- added Sequence::MS_Interface::empty()

version 1.1.2 -- fixed simple bug in output of Sequence::SimParams

version 1.1.1 -- lots of new stuff added to PolyTable.  OS X
compatibility with g++ 3.1 is now working!

version 1.0.9 -- fixed bug in Sequence::PolySites::print() where
position of first seg site was off by 1

version 1.0.8 -- D' added back to
Sequence::Recombination::Disequilibrium()

version 1.0.7 -- removed D' statistic from
Sequence::Recombination::Disequilibrium until a correct version of the code
can be written (I found an error in the calculation and I don't have time to
fix it right now)

version 1.0.6 -- generalized Poly::ThetaW(PolyData) to handle
missing data.  fixed several documentation bugs.

version 1.0.5 -- fixed a problem in SimData.cc that prevented
compiling on Sun systems (did not affect output, just compiling)

version: 1.0.4 -- fixed bug in processing simulation output when
segstites = 0.  Also added a new method to Sequence::SimData that really
speeds up input, see manual for details

version: 1.0.3 -- stable, things work well :)
NOTE: a bug has been found in Sequence::SimData that causes a segfault when
there are 0 segregating sites in the output.  Fix coming soon.



