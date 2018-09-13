# Tutorial/overview

[TOC]

The citation for the library is \cite Thornton2003-wj

## Creation and manipulation of a VariantMatrix

TBW

## Calculation of summary statistics from a VariantMatrix

In libsequence, variation data are represented as a Sequence::VariantMatrix.
The library provides functions for many standard analyses based on input
data in this format.  The following headers are relevant:

1. Sequence/summstats.hpp

Clicking on the above headers will reveal the existence of other headers.
The intent is that you may only wish to bring some names into scope. 
For example, if you implement a new analysis where mean pairwise differences
are needed, you many include Sequence/summstats/thetapi.hpp instead of every 
single summary statistic function provided by the library.

TBW
