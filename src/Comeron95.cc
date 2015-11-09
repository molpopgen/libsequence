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
#include <limits>
#include <algorithm>
#include <Sequence/Seq.hpp>
#include <Sequence/SeqAlphabets.hpp>
#include <Sequence/Comparisons.hpp>
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

/*!
  \defgroup kaks Classes related to the calculation of Ka and Ks
  \ingroup divergence
*/

namespace Sequence
{
  struct Com95impl
  {
    double Qs, Bs, Qa, Ba, A2S, A4, As, A2V, A0, Aa,
      q0, q2S, q2V, q4, p0, p2S, p2V, p4,Ka,Ks;
    std::unique_ptr<RedundancyCom95> sitesObj;
    GeneticCodes code;
    void diverge(const Sequence::Seq & seq1, const Sequence::Seq & seq2,
		 const WeightingScheme2 *_weights2,
		 const WeightingScheme3 *_weights3,
		 const int maxhits);
    void omega (const Sites * s,
		const Sequence::Seq & seqobj1, const Sequence::Seq & seqobj2);
    Com95impl(GeneticCodes __code):
      Qs(0.), Bs(0.), Qa(0.), Ba(0.), A2S(0.), A4(0.), As(0.), A2V(0.), A0(0.), Aa(0.),
      q0(0.), q2S(0.), q2V(0.), q4(0.), p0(0.), p2S(0.), p2V(0.), p4(0.),Ka(0.),Ks(0.),
      sitesObj(std::unique_ptr<RedundancyCom95>(new RedundancyCom95(__code))),
      code(__code)
    {
    }
    
    /*
      assert (max >= 1 && max <=3);
      //check that data are sane
    
      //must be same length
      if(! (seqa.length() == seqb.length()))
	{
	  throw SeqException("Sequence::Comeron95 -- sequences of unequal lengths");
	}
      //must be multiples of 3 in length
      if(!(fabs(double(seqa.length()%3)-0.)<=DBL_EPSILON ) ||
	 !(fabs(double(seqb.length()%3)-0.)<=DBL_EPSILON ) )
	{
	  throw (SeqException("Sequence::Comeron95 -- sequence lengths are not multiples of 3"));
	}

      //If no weighting schemes are passed to the constructor,
      //then we default to using Grantham's distances
      if(_weights2==nullptr)
	{
	  __2wasNULL=true;
	  weights2 = new GranthamWeights2(code);
	}

      if( _weights3 == nullptr)
	{
	  __3wasNULL=true;
	  weights3 = new GranthamWeights3(code);
	}

      maxhits = max;

      q0 = 0.0;
      q2S = 0.0;
      q2V = 0.0;
      q4 = 0.0;
      p0 = 0.0;
      p2S = 0.0;
      p2V = 0.0;
      p4 = 0.0;
    */
    
    double ka (void) const;
    double ks (void) const;
    double ratio (void) const;
    double P0 (void) const;
    double P2S (void) const;
    double P2V (void) const;
    double P4 (void) const;
    double Q0 (void) const;
    double Q2S (void) const;
    double Q2V (void) const;
    double Q4 (void) const;
    double as (void) const;
    double aa (void) const;
    double bs (void) const;
    double ba (void) const;
    double L0 (const Sites *) const;
    double L2S (const Sites * ) const;
    double L2V (const Sites *) const;
    double L4 (const Sites *) const;
  };

  Comeron95::Comeron95( GeneticCodes code ) : impl(std::unique_ptr<Com95impl>(new Com95impl(code)))
  {
  }
									   
    
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

      \exception Sequence::SeqException if sequence lengths are not equal
      \exception Sequence::SeqException if sequence lengths are not multiples of 3
    */

  Com95_t Comeron95::operator()(const Sequence::Seq & seqa,
				const Sequence::Seq & seqb,
				int max)
  {
    GranthamWeights2 w2(impl->code);
    GranthamWeights3 w3(impl->code);
    return this->operator()(seqa,seqb,&w2,&w3,max);
  }
  
  Com95_t Comeron95::operator()(const Sequence::Seq & seqa,
				const Sequence::Seq & seqb,
				const WeightingScheme2 * weights2,
				const WeightingScheme3 * weights3,
				int max)
  {
    Sites s(impl->sitesObj.get(),&seqa,&seqb,max,impl->code);
    impl->diverge(seqa,seqb,weights2,weights3,max);
    impl->omega(&s,seqa,seqb);
    return Com95_t({impl->ka(),
	  impl->ks(),
	  impl->ratio(),
	  impl->P0(),
	  impl->P2S(),
	  impl->P2V(),
	  impl->P4(),
	  impl->Q0(),
	  impl->Q2S(),
	  impl->Q2V(),
	  impl->Q4(),
	  impl->as(),
	  impl->aa(),
	  impl->bs(),
	  impl->ba(),
	  impl->L0(&s),
	  impl->L2S(&s),
	  impl->L2V(&s),
	  impl->L4(&s)
	  });
  }
  
  void Com95impl::omega (const Sites * s,
			 const Sequence::Seq & seqobj1,
                         const Sequence::Seq & seqobj2)
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

    Qs = (q2V + q4) / (s->L2V() + s->L4());

    if (!std::isfinite (Qs))
      Qs = 0.0;

    Bs = (-0.5) * log (1.0 - (2.0 * Qs));

    Qa = (q0 + q2S) / (s->L0() + s->L2S());

    if (!std::isfinite (Qa))
      Qa = 0.0;

    Ba = (-0.5) * std::log (1.0 - (2.0 * Qa));
    if (!std::isfinite (Ba))// && !isinf(Ba))
      {
        //if sites aren't saturated in general (i.e. Kimura's distance < 1.0)
        //it is likely that Ba is nan due to too few changes, and thus Ba should equal 0.0
        //otherwise it is due to too many changes.
        //NOTE--this is an ad-hoc treatment of the analysis!!!
        std::unique_ptr<Kimura80> K80( new Kimura80 (&seqobj1, &seqobj2));
        if (K80->K () < 1.0)
          Ba = 0.0;
      }
    
    //calculate numbers of mutation per site type
    double P2S_site = p2S / s->L2S();
    double P2V_site = p2V / s->L2V();
    double P0_site = p0 / s->L0();
    double Q0_site = q0 / s->L0();
    double P4_site = p4 / s->L4();
    double Q4_site = q4 / s->L4();

    log1 = std::log (1.0 - (2.0 * P2S_site) - Qa);
    log2 = std::log (1.0 - (2.0 * Qa));

    if (!std::isfinite (log1))	//set value to 0 if nan
      log1 = 0.0;
    if (!std::isfinite (log2))
      log2 = 0.0;

    A2S = (-0.5) * log1 + (0.25) * log2;

    log1 = std::log (1.0 - (2.0 * P4_site) - Q4_site);
    log2 = std::log (1.0 - (2.0 * Q4_site));

    if (!std::isfinite (log1))
      log1 = 0.0;
    if (!std::isfinite (log2))
      log2 = 0.0;

    A4 = (-0.5) * log1 + (0.25) * log2;

    As = (s->L2S() * A2S + s->L4() * A4) / (s->L2S() +
						    s->L4());

    log1 = std::log (1.0 - (2.0 * P2V_site) - Qs);
    log2 = std::log (1.0 - (2.0 * Qs));

    if (!std::isfinite (log1))
      log1 = 0.0;
    if (!std::isfinite (log2))
      log2 = 0.0;

    A2V = (-0.5) * log1 + (0.25) * log2;

    log1 = std::log (1.0 - (2.0 * P0_site) - Q0_site);
    log2 = std::log (1.0 - (2.0 * Q0_site));
    if (!std::isfinite (log1))
      log1 = 0.;
    if (!std::isfinite (log2))
      log2 = 0.0;

    A0 = (-0.5) * log1 + (0.25) * log2;

    Aa = (s->L2V() * A2V + s->L0() * A0) / (s->L2V() +
						    s->L0());

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

    if (!std::isfinite (Ks))
      Ks = std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite (Ka))
      Ka = std::numeric_limits<double>::quiet_NaN();
  }

  void Com95impl::diverge (const Sequence::Seq & seq1,
                           const Sequence::Seq & seq2,
                           const WeightingScheme2 *weights2,
                           const WeightingScheme3 *weights3,
			   const int maxhits)
  /*!
    go through every aligned, ungapped codon,
    and calculate divergence.  maintains a running sum of divergence
    statistics stored a private data to the class
  */
  {
    size_t i, j;
    size_t length = seq1.length();
     q0= q2S= q2V= q4= p0= p2S= p2V= p4 = 0.;
    //the for loop iterates over codons (block of 3 sites)
    std::string codon1, codon2;
    codon1.resize (3);
    codon2.resize (3);

    for (i = 0; i <= (length - 3); i += 3)
      {
        for (j = 0; j <= 2; ++j)
          {
            //assign the next codon from each sequence
            codon1[j] = char(std::toupper(seq1[i + j]));
            codon2[j] = char(std::toupper(seq2[i + j]));
          }
        if (  std::find_if(codon1.begin(),codon1.end(),ambiguousNucleotide()) == codon1.end() &&
	      std::find_if(codon2.begin(),codon2.end(),ambiguousNucleotide()) == codon2.end() )
          {
            //find out if codons are different
            if (Different (codon1, codon2))
              {
                //find out how many changes there are between the codons
                int ndiff = NumDiffs (codon1, codon2);
                if (ndiff == 1)	//if there is just one difference, use the rules for
                  //codons that only differ at 1 site
                  {
                    SingleSub Single;
                    Single(sitesObj.get(), codon1, codon2);
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
                            Double(sitesObj.get(), codon1, codon2,weights2);
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
                            Triple(sitesObj.get(),codon1, codon2,weights3);
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

    if (!std::isfinite (p0))
      p0 = 0.0;
    if (!std::isfinite (p2S) )
      p2S = 0.0;
    if (!std::isfinite (p2V) )
      p2V = 0.0;
    if (!std::isfinite (p4) )
      p4 = 0.0;
    if (!std::isfinite (q0) )
      q0 = 0.0;
    if (!std::isfinite (q2S) )
      q2S = 0.0;
    if (!std::isfinite (q2V) )
      q2V = 0.0;
    if (!std::isfinite (q4) )
      q4 = 0.0;
  }

  double Com95impl::L0 (const Sites * sites) const
  /*!
    \return the number of nondegenerate sites compared
  */
  {
    return sites->L0();
  }
  double Com95impl::L2S (const Sites * sites) const
  /*!
    \return the number of twofold, transitional-degenerate sites compared
  */
  {
    return sites->L2S();
  }
  double Com95impl::L2V (const Sites * sites) const
  /*!
    \return the number of twofold, transversional-degenerate sites compared
  */
  {
    return sites->L2V();
  }
  double Com95impl::L4 (const Sites * sites) const
  /*!
    \return the number of 4-fold degenerate sites compared
  */
  {
    return sites->L4();
  }

  double Com95impl::as (void) const
  /*!
    \return corrected synonymous divergence at transitional-degenerate sites
  */
  {
    if (!std::isfinite ( As))
      return std::numeric_limits<double>::quiet_NaN();
    return As;
  }
  double Com95impl::aa (void) const
  /*!
    \return corrected nonsynonymous divergence at tranversioal- and non- degenerate sites
  */
  {
    if (!std::isfinite ( Aa))
      return std::numeric_limits<double>::quiet_NaN();
    return Aa;
  }
  double Com95impl::bs (void) const
  /*!
    \return corrected synonymous divergence at transversional- and fourfold-  degenerate sites
  */
  {
    if (!std::isfinite ( Bs))
      return std::numeric_limits<double>::quiet_NaN();
    return Bs;
  }
  double Com95impl::ba (void) const
  /*!
    \return corrected nonsynonymous divergence at transitional- and non- degenerate sites
  */
  {
    if (!std::isfinite ( Ba))
      return std::numeric_limits<double>::quiet_NaN();
    return Ba;
  }

  double Com95impl::ratio(void) const
  /*!
    \return \f$K_a/K_s\f$
    \note std::numeric_limits<double>::quiet_NaN() is returned if Ka/Ks cannot be calculated
  */
  {
    if ( !std::isfinite(Ka) || !std::isfinite(Ks) || fabs(Ks-0.) <= DBL_EPSILON)
      return std::numeric_limits<double>::quiet_NaN();
    return Ka / Ks;
  }
  
  double Com95impl::ka (void) const
  /*!
    \return the nonsynonymous distance
    \note std::numeric_limits<double>::quiet_NaN() is returned if Ka cannot be calculated
  */
  {
    return Ka;
  }
  double Com95impl::ks (void) const
  /*!
    \return the synonymous distance
    \note std::numeric_limits<double>::quiet_NaN() is returned if Ks cannot be calculated
  */
  {
    return Ks;
  }

  double Com95impl::P0 (void) const
  /*!
    \return number of transitions at nondegenerate sites
  */
  {
    return p0;
  }
  double Com95impl::P2S (void) const
  /*!
    \return number of transitions at 2-fold, transitional degenerate sites
  */
  {
    return p2S;
  }
  double Com95impl::P2V (void) const
  /*!
    \return number of transitions at  2-fold, transversional degenerate sites
  */
  {
    return p2V;
  }
  double Com95impl::P4 (void) const
  /*!
    \return number of transitions at 4-fold degenerate sites
  */
  {
    return p4;
  }
  double Com95impl::Q0 (void) const
  /*!
    \return number of transversion at nondegenerate sites
  */
  {
    return q0;
  }
  double Com95impl::Q2S (void) const
  /*!
    \return number of transversion at 2-fold, transitional degenerate sites
  */
  {
    return q2S;
  }
  double Com95impl::Q2V (void) const
  /*!
    \return number of transversion at 2-fold, transversional sites
  */
  {
    return q2V;
  }
  double Com95impl::Q4 (void) const
  /*!
    \return number of transversion at 4-fold degenerate sites
  */
  {
    return q4;
  }
}
