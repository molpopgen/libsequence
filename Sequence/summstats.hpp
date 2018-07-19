/// @file Sequence/summstats.hpp
/// \brief Include all summary statistic functions and types
#ifndef SEQUENCE_SUMMSTATS_HPP__
#define SEQUENCE_SUMMSTATS_HPP__

/*!
 *  \defgroup popgenanalysis Analysis of molecular population genetic data
 *  \brief Summary statistics and other analysis of Sequence::VariantMatrix
 *  \ingroup popgen

 * In libsequence, variation data are represented as a Sequence::VariantMatrix.
 * The library provides functions for many standard analyses based on input
 * data in this format.  The following headers are relevant:

 * 1. Sequence/summstats.hpp
 *
 * Clicking on the above headers will reveal the existence of other headers.
 * The intent is that you may only wish to bring some names into scope. 
 * For example, if you implement a new analysis where mean pairwise differences
 * are needed, you many include Sequence/summstats/thetapi.hpp instead of every 
 * single summary statistic function provided by the library.
 * 
*/

#include "summstats/classics.hpp"

#endif
