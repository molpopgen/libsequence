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

#ifndef RECOMBINATION_H
#define RECOMBINATION_H

/*! \file Sequence/Recombination.hpp
  @brief namespace Sequence::Recombination
*/

/*! \namespace Sequence::Recombination
  \ingroup popgenanalysis
  This namespace exists primarily so that the file Poly.cc (which defines
  Sequence::Poly) does not get too large.  The routines defined in this
  namespace
  all have to do with properties of the association between sites.  Current
  methods implemented
  are:\n
  1.) Recombination::HudsonsC, which calculates Hudson's C, aka
  \f$\rho_{87}\f$\n
  \n
  2.) Recombination::Disequilibrium, which calculated several measures of LD
  for all pairs of sites,
  and implements a frequency filter to remove low-frequency variants if
  desired.

  @short Methods dealing with recombination
*/
#include <vector>
#include <limits>

namespace Sequence
{
    struct PairwiseLDstats
		/*! \brief Pairwise linkage disequilibrium (LD) stats
		 * \ingroup popgenanalysis
		 */
    {
		//! Position of site i
        double i;
		//! Position of site j
		double j;
		//! $r^2$ between sites i and j
		double rsq;
		//! The D statistics
		double D;
		//! D'
		double Dprime;
		//! If the site fails and filters, this is true
        bool skipped;
        PairwiseLDstats();
    };

    class PolyTable;
    namespace Recombination
    {
        double HudsonsC(const Sequence::PolyTable *data,
                        const bool &haveOutgroup, const unsigned &outgroup);
        /*!
		 * \brief Calculate pairwise LD for a Sequence::PolyTable
		 * \param data A Sequence::PolyTable
		 * \param haveOutgroup A boolean
		 * \param outgroup. If \a haveOutgroup is true, then this is
         * the index
		 * of the outgroup sequence in \a data
		 * \param mincount  Minimum sample count to include a mutation.
         * E.g., 2 =
		 * exclude singletons, etc.
		 * \param max_distance Do not include sites > this value apart
		 * \return A std::vector of Sequence::PairwiseLDstats
		 */
        std::vector<PairwiseLDstats> Disequilibrium(
            const Sequence::PolyTable *data, const bool &haveOutgroup = false,
            const unsigned &outgroup = 0, const unsigned &mincount = 1,
            const double max_distance = std::numeric_limits<double>::max());

        /*!
         Calculates LD statistics for sites i and j, where j>1;
		 \param data The polymorphism data
         \param i Index for site in data
         \param j Index for site in data
         \param haveOutgroup true if \a data contains an outgroup, false
         otherwise
         \param outgroup the index of the outgroup in \a data \a haveOutgroup
         is true
         \param mincount the minimum sample frequency to include
         \param max_distance max distance between markers (on scale in which
         positions in \a data are stored) for inclusion
		 \return Sequence::PairwiseLDstats
       */
        PairwiseLDstats PairwiseLD(const Sequence::PolyTable *data, unsigned i,
                                   unsigned j,
                                   const bool &haveOutgroup = false,
                                   const unsigned &outgroup = 0,
                                   const unsigned &mincount = 1,
                                   const double max_distance
                                   = std::numeric_limits<double>::max());
    }
}
#endif
