#What this document implies

This document lists major changes between tags/releases of libsequence. The following caveat is important:

If you see a section below referring to libseqeunce version X, and no version X exists as a release, then you must consider the info here inaccurate.  In such cases, the list of changes may refer to planned changes for a future release which may or may not be implemented yet in either the master or some development branch.

##libsequence 1.8.3

(1.8.3 is unreleased.  The list below may be incomplete, suffer from wishful thinking, etc.)

1.  This release _breaks binary compatibility_ with previous releases.  Programs depending on libsequence will likely need to be recompiled.
2.  New constructor functions were added to Sequence::PolyTable, which breaks binary compatibility.  These new constructors are needed for easier [Rcpp](http://www.rcpp.org/) integration.  Derived classes have been updated to use these new constructors whenever possible.
3.  Some I/O functions are no longer inline.  These include the operator<< and operator>> for Sequence::Seq and Sequence::PolyTable.
