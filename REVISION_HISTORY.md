#What this document implies

This document lists major changes between tags/releases of libsequence. The following caveat is important:

If you see a section below referring to libseqeunce version X, and no version X exists as a release, then you must consider the info here inaccurate.  In such cases, the list of changes may refer to planned changes for a future release which may or may not be implemented yet in either the master or some development branch.

##libsequence 1.8.3

(1.8.3 is unreleased.  The list below may be incomplete, suffer from wishful thinking, etc.)

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

###ISSUES:

1. Sequence/SeqRegexes.hpp -- not working.  This will not be fixed until GCC supports <regex>.  The function is now currently implemented in a non-regex manner, which is lame, but it works.

###TODO:

1. H12 from Petrov, Messer, et al.
2. G stat for differentiation -- http://arxiv.org/pdf/1403.1552.pdf
3. iHH?
4. http://mbe.oxfordjournals.org/content/31/5/1275.abstract