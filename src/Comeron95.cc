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

#include <cmath>
#include <cfloat>
#include <cctype>
#include <memory>
#include <cassert>
#include <algorithm>
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#include <Sequence/Seq.hpp>
#include <Sequence/SeqProperties.hpp>
#include <Sequence/Grantham.hpp>
#include <Sequence/GranthamWeights.hpp>
#include <Sequence/SingleSub.hpp>
#include <Sequence/TwoSubs.hpp>
#include <Sequence/ThreeSubs.hpp>
#include <Sequence/Sites.hpp>
#include <Sequence/Kimura80.hpp>
#include <Sequence/RedundancyCom95.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <Sequence/Comeron95.hpp>

using std::log;
using std::toupper;
using std::isfinite;
/*!
  \defgroup kaks Classes related to the calculation of Ka and Ks
  \ingroup divergence
*/

namespace Sequence
  {
  Comeron95::Comeron95 (const Sequence::Seq * seqa,
                        const Sequence::Seq * seqb, 
			int max, 
			const Sequence::RedundancyCom95 * genetic_code_redundancy,
			GeneticCodes code,
                        WeightingScheme2 *_weights2,
                        WeightingScheme3 *_weights3) :
    __2wasNULL(false),__3wasNULL(false),__red_was_NULL(false)
  /*!
    Initialize and calculate synonymous and nonsynonymous distances between two sequence objects
    \param seqa an object of type or derived from type Sequence::Seq
    \param seqb an object of type or derived from type Sequence::Seq
    \param max maximum number of substitutions per codon to allow in the analysis
    \param code genetic code, see Sequence::GeneticCodes
    \param _weights2 a weighting scheme for codons differing at 2 positions.  
    If \c NULL, Sequence::GranthamWeights2 is used
    \param _weights3 a weighting scheme for codons differing at 3 positions.  
    If \c NULL, Sequence::GranthamWeights3 is used
    \warning Note that the pointers to weighting schemes are dumb pointers.
    This allows me to check for NULL and then assign a default.  If you
    use your own classes, make sure they clean up after themselves if they
    throw exceptions!!!
    \todo use \c auto_ptr for the weighting schemes
    \exception Sequence::SeqException if sequence lengths are not equal
    \exception Sequence::SeqException if sequence lengths are not multiples of 3
  */
  {
    assert (max >= 1 && max <=3);
    //check that data are sane
    
    //must be same length
    if(! (seqa->length() == seqb->length()))
      {
        throw SeqException("Sequence::Comeron95 -- sequences of unequal lengths");
      }
    //must be multiples of 3 in length
    if(!(fabs(double(seqa->length()%3)-0.)<=DBL_EPSILON ) ||
       !(fabs(double(seqb->length()%3)-0.)<=DBL_EPSILON ) )
      {
        throw (SeqException("Sequence::Comeron95 -- sequence lengths are not multiples of 3"));
      }

    //If no weighting schemes are passed to the constructor,
    //then we default to using Grantham's distances
    if(_weights2==NULL)
      {
        __2wasNULL=true;
        weights2 = new GranthamWeights2(code);
      }

    if( _weights3 == NULL)
      {
        __3wasNULL=true;
        weights3 = new GranthamWeights3(code);
      }

    maxhits = max;
    if (genetic_code_redundancy == NULL)
      {
	sitesObj = new RedundancyCom95 (code);
	__red_was_NULL = true;
      }
    else
      sitesObj = genetic_code_redundancy;

    sites = new Sites (sitesObj, seqa, seqb, maxhits, code);
    q0 = 0.0;
    q2S = 0.0;
    q2V = 0.0;
    q4 = 0.0;
    p0 = 0.0;
    p2S = 0.0;
    p2V = 0.0;
    p4 = 0.0;

    diverge (seqa, seqb,weights2,weights3);
    omega (seqa, seqb);
  }

  Comeron95::~Comeron95 (void)
  /*!
    deletes a pointer to a Sequence::Sites and a Sequence::RedundancyCom95
  */
  {
    delete sites;
    //if we had to initialize weighting schemes,
    //we need to delete them here to
    //prevent memory leaks...
    if(__2wasNULL==true)
      delete weights2;
    if(__3wasNULL==true)
      delete weights3;
    if (__red_was_NULL==true)
      delete sitesObj;
  }

  void Comeron95::omega (const Sequence::Seq * seqobj1,
                         const Sequence::Seq * seqobj2)
  /*!
    calculate values needed to obtain Ka and Ks.
    formulae are from Comeron '95 and use the identical notation
    the rest of the code is just calculating numbers from the data.
    At first glance, it looks like this is a lot of code,
    but these calculations require careful checking, becuase when you
    have a lot of changes, or very few, you can take logs of negative
    numbers, or divide by zero, hence calls to isnan() and isinf()
  */
  {
    double log1, log2;

    Qs = (q2V + q4) / (sites->L2V() + sites->L4());

    if (!isfinite (Qs))
      Qs = 0.0;

    Bs = (-0.5) * log (1.0 - (2.0 * Qs));

    //     if (isnan (Bs))
    //       {
    // 	Kimura80 *K80 = new Kimura80 (seqobj1, seqobj2);
    // 	//if sites aren't saturated in general (i.e. Kimura's distance < 1.0)
    // 	//it is likely that Bs is nan due to too few changes, and thus Bs should equal 0.0
    // 	//otherwise it is due to too many changes.
    // 	//NOTE--this is an ad-hoc treatment of the analysis!!!
    // 	if (K80->K () < 1.0)
    // 	  Bs = 0.0;
    // 	delete (K80);
    //       }

    Qa = (q0 + q2S) / (sites->L0() + sites->L2S());

    if (!isfinite (Qa))
      Qa = 0.0;

    Ba = (-0.5) * log (1.0 - (2.0 * Qa));
    if (!isfinite (Ba))// && !isinf(Ba))
      {
        //if sites aren't saturated in general (i.e. Kimura's distance < 1.0)
        //it is likely that Ba is nan due to too few changes, and thus Ba should equal 0.0
        //otherwise it is due to too many changes.
        //NOTE--this is an ad-hoc treatment of the analysis!!!
        std::auto_ptr<Kimura80> K80( new Kimura80 (seqobj1, seqobj2));
        if (K80->K () < 1.0)
          Ba = 0.0;
      }

    //calculate numbers of mutation per site type
    double P2S_site = p2S / sites->L2S();
    double P2V_site = p2V / sites->L2V();
    double P0_site = p0 / sites->L0();
    double Q0_site = q0 / sites->L0();
    double P4_site = p4 / sites->L4();
    double Q4_site = q4 / sites->L4();

    log1 = log (1.0 - (2.0 * P2S_site) - Qa);
    log2 = log (1.0 - (2.0 * Qa));

    if (!isfinite (log1))	//set value to 0 if nan
      log1 = 0.0;
    if (!isfinite (log2))
      log2 = 0.0;

    A2S = (-0.5) * log1 + (0.25) * log2;

    log1 = log (1.0 - (2.0 * P4_site) - Q4_site);
    log2 = log (1.0 - (2.0 * Q4_site));

    if (!isfinite (log1))
      log1 = 0.0;
    if (!isfinite (log2))
      log2 = 0.0;

    A4 = (-0.5) * log1 + (0.25) * log2;

    As = (sites->L2S() * A2S + sites->L4() * A4) / (sites->L2S() +
         sites->L4());

    log1 = log (1.0 - (2.0 * P2V_site) - Qs);
    log2 = log (1.0 - (2.0 * Qs));

    if (!isfinite (log1))
      log1 = 0.0;
    if (!isfinite (log2))
      log2 = 0.0;

    A2V = (-0.5) * log1 + (0.25) * log2;

    log1 = log (1.0 - (2.0 * P0_site) - Q0_site);
    log2 = log (1.0 - (2.0 * Q0_site));
    if (!isfinite (log1))
      log1 = 0.;
    if (!isfinite (log2))
      log2 = 0.0;

    A0 = (-0.5) * log1 + (0.25) * log2;

    Aa = (sites->L2V() * A2V + sites->L0() * A0) / (sites->L2V() +
         sites->L0());

    if (As <= 0.0)
      As = 0.0;
    if (Bs <= 0.0)
      Bs = 0.0;
    if (Aa <= 0.0)
      Aa = 0.0;
    if (Ba <= 0.0)
      Ba = 0.0;

    Ks = As + Bs;
    Ka = Aa + Ba;

    if (!isfinite (Ks))
      Ks = 999;
    if (!isfinite (Ka))
      Ka = 999;
  }

  void Comeron95::diverge (const Sequence::Seq * seq1,
                           const Sequence::Seq * seq2,
                           WeightingScheme2 *weights2,
                           WeightingScheme3 *weights3)
  /*!
    go through every aligned, ungapped codon,
    and calculate divergence.  maintains a running sum of divergence
    statistics stored a private data to the class
  */
  {
    size_t i, j, ndiff;
    size_t length = seq1->length ();

    //the for loop iterates over codons (block of 3 sites)
    std::string codon1, codon2;
    codon1.resize (3);
    codon2.resize (3);

    for (i = 0; i <= (length - 3); i += 3)
      {
        for (j = 0; j <= 2; ++j)
          {
            //assign the next codon from each sequence
            codon1[j] = char(std::toupper((*seq1)[i + j]));
            codon2[j] = char(std::toupper((*seq2)[i + j]));
          }
        if (  std::find_if(codon1.begin(),codon1.end(),ambiguousNucleotide()) == codon1.end() &&
	      std::find_if(codon2.begin(),codon2.end(),ambiguousNucleotide()) == codon2.end() )
          {
            //find out if codons are different
            if (Different (codon1, codon2))
              {
                //find out how many changes there are between the codons
                ndiff = NumDiffs (codon1, codon2);
                if (ndiff == 1)	//if there is just one difference, use the rules for
                  //codons that only differ at 1 site
                  {
                    SingleSub Single;
                    Single(sitesObj, codon1, codon2);
                    p0 += Single.P0 ();
                    p2S += Single.P2S ();
                    p2V += Single.P2V ();
                    p4 += Single.P4 ();
                    q0 += Single.Q0 ();
                    q2S += Single.Q2S ();
                    q2V += Single.Q2V ();
                    q4 += Single.Q4 ();
                  }

                if (maxhits > 1)
                  {	//if codons with >1 difference are allowed
                    //iterate over codons, as above

                    if (Different
                        (codon1, codon2))
                      {
                        ndiff = NumDiffs
                                (codon1,
                                 codon2);
                        //cout << "CHECK:NDIFF="<<ndiff<<endl;
                        if (ndiff == 2
                            && maxhits >= 2)
                          {	//if codons differ at 2 sites,
                            //use rules of class TwoSubstitutions
                            TwoSubs Double;
                            Double(sitesObj, codon1, codon2,weights2);
                            p0 += Double. P0();
                            p2S += Double.P2S ();
                            p2V += Double.P2V ();
                            p4 += Double.P4();
                            q0 += Double.Q0 ();
                            q2S += Double.Q2S ();
                            q2V += Double.Q2V ();
                            q4 += Double.Q4 ();
                          }
                        else if (ndiff == 3 && maxhits > 2)
                          {
                            ThreeSubs Triple;
                            Triple(sitesObj,codon1, codon2,weights3);
                            p0 += Triple. P0 ();
                            p2S += Triple.P2S ();
                            p2V += Triple.P2V ();
                            p4 += Triple. P4 ();
                            q0 += Triple.Q0 ();
                            q2S += Triple.Q2S ();
                            q2V += Triple.Q2V ();
                            q4 += Triple.Q4 ();
                          }
                      }
                  }
              }
          }
      }

    if (!isfinite (p0))
      p0 = 0.0;
    if (!isfinite (p2S) )
      p2S = 0.0;
    if (!isfinite (p2V) )
      p2V = 0.0;
    if (!isfinite (p4) )
      p4 = 0.0;
    if (!isfinite (q0) )
      q0 = 0.0;
    if (!isfinite (q2S) )
      q2S = 0.0;
    if (!isfinite (q2V) )
      q2V = 0.0;
    if (!isfinite (q4) )
      q4 = 0.0;
  }

  double Comeron95::L0 (void) const
  /*!
    \return the number of nondegenerate sites compared
  */
  {
    return sites->L0();
  }
  double Comeron95::L2S (void) const
  /*!
    \return the number of twofold, transitional-degenerate sites compared
  */
  {
    return sites->L2S();
  }
  double Comeron95::L2V (void) const
  /*!
    \return the number of twofold, transversional-degenerate sites compared
  */
  {
    return sites->L2V();
  }
  double Comeron95::L4 (void) const
  /*!
    \return the number of 4-fold degenerate sites compared
  */
  {
    return sites->L4();
  }

  double Comeron95::as (void) const
  /*!
    \return corrected synonymous divergence at transitional-degenerate sites
  */
  {
    if (!isfinite (As))
      return 999.;
    return As;
  }
  double Comeron95::aa (void) const
  /*!
    \return corrected nonsynonymous divergence at tranversioal- and non- degenerate sites
  */
  {
    if (!isfinite (Aa))
      return 999.;
    return Aa;
  }
  double Comeron95::bs (void) const
  /*!
    \return corrected synonymous divergence at transversional- and fourfold-  degenerate sites
  */
  {
    if (!isfinite (Bs))
      return 999.;
    return Bs;
  }
  double Comeron95::ba (void) const
  /*!
    \return corrected nonsynonymous divergence at transitional- and non- degenerate sites
  */
  {
    if (!isfinite (Ba))
      return 999.;
    return Ba;
  }

  double Comeron95::ratio(void) const
  /*!
    \return \f$K_a/K_s\f$
    \note 999.0 is returned if Ka/Ks cannot be calculated
  */
  {

    if (Ka == 999. || Ks == 999. || fabs(Ks-0.) <= DBL_EPSILON)
      return 999.;
    return Ka / Ks;
  }
  double Comeron95::ka (void) const
  /*!
    \return the nonsynonymous distance
    \note 999.0 is returned if Ka cannot be calculated
  */
  {
    return Ka;
  }
  double Comeron95::ks (void) const
  /*!
    \return the synonymous distance
    \note 999.0 is returned if Ks cannot be calculated
  */
  {
    return Ks;
  }

  double Comeron95::P0 (void) const
  /*!
    \return number of transitions at nondegenerate sites
  */
  {
    return p0;
  }
  double Comeron95::P2S (void) const
  /*!
    \return number of transitions at 2-fold, transitional degenerate sites
  */
  {
    return p2S;
  }
  double Comeron95::P2V (void) const
  /*!
    \return number of transitions at  2-fold, transversional degenerate sites
  */
  {
    return p2V;
  }
  double Comeron95::P4 (void) const
  /*!
    \return number of transitions at 4-fold degenerate sites
  */
  {
    return p4;
  }
  double Comeron95::Q0 (void) const
  /*!
    \return number of transversion at nondegenerate sites
  */
  {
    return q0;
  }
  double Comeron95::Q2S (void) const
  /*!
    \return number of transversion at 2-fold, transitional degenerate sites
  */
  {
    return q2S;
  }
  double Comeron95::Q2V (void) const
  /*!
    \return number of transversion at 2-fold, transversional sites
  */
  {
    return q2V;
  }
  double Comeron95::Q4 (void) const
  /*!
    \return number of transversion at 4-fold degenerate sites
  */
  {
    return q4;
  }
}
