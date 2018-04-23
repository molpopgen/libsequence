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

#ifndef _POLYSNP_H_
#define _POLYSNP_H_
/*! \file PolySNP.hpp
  @brief declaration of Sequence::PolySNP, a class to analyze SNP data
*/
/*!
  \class Sequence::PolySNP Sequence/PolySNP.hpp
  \ingroup popgenanalysis
  Example Use:\n
  \n
  \code
  #include <iostream>
  #include <vector>
  #include <Sequence/PolySNP.hpp>
  #include <Sequence/Fasta.hpp>
  #include <Sequence/PolySites.hpp>
  using namespace std;
  using namespace Sequence;
  int main
  {
    vector<Fasta> data;
    Alignment::GetData(data,"popdata.fasta");
    assert(Alignment::IsAlignment(data))
    if (Alignment::Gapped(data))
    {
      Alignment::RemoveTerminalGaps(data);
    }
    PolySites *polytable = new PolySites(data);
    PolySNP *analyze = new PolySNP(data,false,0);
    cout << "Tajima's D is " << analyze->TajimasD() << endl;
    delete polytable;
    delete analyze;
    exit(1);
  }
  \endcode
  \warning The routines to calculate nucleotide diversity, numbers of singletons, etc. all handle
  missing data (the \c N character).  However, all summary statistics involved in "tests of 
  neutrality" are, strictly speaking, undefined if missing data are present.  The reason for this is
  because the denominators of the statistics are functions of the sample sizes, and no explicit
  formulae exist when the sample size varies from site to site (which is the case when there are
  missing data).  In short, if you want to be rigorous, you can only really count up nucleotide diversity
  and a few other statistics if your data contain untyped SNPs.  However, the routines present in
  \c libsequence will happily go and calculate the summary statistics for you, and it is up to you
  to be aware that you are writing a program that may give biased results.  To date, the magnitude
  and direction of the bias remains unknown.  Functions (and hence the statistics) that 
  are affected have warnings in their documentation.
  \note As of libsequence 1.4.1, the routines in this class explicity check for gaps in the
  polymorphism table.  This provides an additional safeguard for cases where Sequence::PolyTable
  objects are created and sites with gaps are left in.
  @short Molecular population genetic analysis

  \deprecated Will be removed in libsequence 2.0
*/
#include <vector>
#include <memory>
#include <limits>
#include <Sequence/Recombination.hpp>

namespace Sequence 
  {
    class __attribute__ ((deprecated)) PolyTable;
    struct _PolySNPImpl;
    class PolySNP 
    {
    private:
    protected:
      std::unique_ptr<_PolySNPImpl> rep;
      void DepaulisVeuilleStatistics (void) const;
      virtual void WallStats(void) const;
      //various things one needs to know to calculate the summary statistics
      double a_sub_n (void) const ;
      double a_sub_n_plus1 (void) const;
      double b_sub_n (void) const;
      double b_sub_n_plus1 (void) const;
      double c_sub_n (void) const;
      double d_sub_n (void) const;
    public:
      explicit PolySNP (const Sequence::PolyTable * data, bool haveOutgroup = false,
                        unsigned outgroup = 0, bool totMuts = true);
      PolySNP( const PolySNP & ) = delete;
      PolySNP & operator=(const PolySNP &) = delete;
      virtual ~ PolySNP (void);
      //estimators of 4Nu
      virtual double ThetaPi (void) const;                           //Nucleotide diversity (Tajima 1983)
      virtual double ThetaW (void) const;                            //Watterson's (1975) Theta
      virtual double ThetaH (void) const;                            //Theta from homozygosity (Fay and Wu 2001)
      virtual double ThetaL (void) const;                            //A variant on Fay and Wu's H
      //variances of estimators of 4Nu
      double VarPi (void) const;
      double StochasticVarPi(void) const;
      double SamplingVarPi (void) const;
      double VarThetaW (void) const;
      //calculate various numbers related to polymorphism
      unsigned NumPoly (void) const;                                 //number of polymorphic sites in data
      virtual unsigned NumMutations (void) const;                    //number of inferred mutations 
      virtual unsigned NumSingletons (void) const;                   //number of mutants at frequency 1
      virtual unsigned NumExternalMutations (void) const;            //number of derived mutations 
      //summary statistics of the site frequency spectrum
      virtual double TajimasD (void) const;                          //Tajima's (1989) D
      virtual double Hprime (const bool & likeThorntonAndolfatto = false) const;  //A normalized statistic 
                                                                    //related to Fay and Wu's H
      virtual double Dnominator (void) const;                        //Denominator of Tajima's D
      virtual double FuLiD (void) const;                             //Fu & Li's (1996) D
      virtual double FuLiF (void) const;                             //Fu & Li's (1996) F
      virtual double FuLiDStar (void) const;                         //Fu & Li's (1996) D*
      virtual double FuLiFStar (void) const;                         //Fu & Li's (1996) F*
      //summary statistics of haplotypes
      double DandVH (void) const;                                    //Depaulis & Veuille (1998) Haplotype diversity
      unsigned DandVK (void) const;                                  //Depaulis & Veuille (1998) number of haplotypes
      virtual double WallsB(void) const;                             //Jeff Wall's (1999) B  statistic
      virtual unsigned WallsBprime(void) const;                      //Jeff Wall's (1999) B' statistic
      virtual double WallsQ(void) const;                             //Jeff Wall's (1999) Q  statistic
      //recombination
      double HudsonsC (void) const;                                  //Dick Hudson's (1987) Chat = 4Nr
      virtual unsigned Minrec (void) const;                          //Hudson & Kaplan's (1985) min # of recombination events
      std::vector <PairwiseLDstats>
      Disequilibrium( const unsigned & mincount = 1,
		      const double & max_marker_distance = std::numeric_limits<double>::max()) const;                 //summary stats of pairwise LD
    };
}
#endif
