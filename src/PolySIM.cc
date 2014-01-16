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
#include <cassert>
#include <cstdlib>
#include <Sequence/SimData.hpp>
#include <Sequence/Recombination.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolySNPimpl.hpp>
/*!
  \ingroup popgenanalysis
*/

using namespace Sequence::Recombination;

namespace Sequence
{
  PolySIM::PolySIM (const Sequence::SimData * data):
    PolySNP(data,false,0,true)
		   /*!
		     \param data a valid object of type Sequence::SimData
		   */
  {
    rep->_NumPoly = NumPoly();
  }

  PolySIM::~PolySIM(void)
  {
  }

  double
  PolySIM::ThetaPi (void)
  /*!
    For simulated data, assuming 0 is ancenstral, 1 derived.\n
    A simpler version of PolySNP::ThetaPi
  */
  {
    if (rep->_know_pi == false)
      {
        double Pi = 0.0, nsam = double(rep->_nsam);
	for( std::vector<stateCounter>::const_iterator i = rep->_counts.begin() ;
	     i != rep->_counts.end() ; ++i)
	  {
	    Pi += (2.0 * i->one *(nsam-i->one)) / (nsam*(nsam-1.));
	  }
        rep->_pi = Pi;
        rep->_know_pi = true;
      }
    return rep->_pi;
  }

  double
  PolySIM::ThetaW (void)
  /*!
    For coalescent simulation data, the number of segregating 
    sites equals the number of mutations
    on the tree (under the infinite sites model.  
  */
  {
    return (rep->_NumPoly>0) ? (double(rep->_NumPoly)/a_sub_n() ) : 0.;
  }


  double
  PolySIM::ThetaH (void)
  /*!
    For simulated data, where 0 is ancenstral, 1 derived.\n
    A simpler version of PolySIM::ThetaH (const Sequence::PolyTable * data,
    bool haveOutgroup = 0, unsigned outgroup = 0)
  */
  {
    double H = 0.0, nsam = double(rep->_nsam);

    for( std::vector<stateCounter>::const_iterator i = rep->_counts.begin();
	 i != rep->_counts.end() ; ++i)
      {
	H += (i->one < nsam) ? (2.0 * i->one * i->one)/ (nsam * (nsam - 1.0)) : 0.;
      }
    return H;
  }

  double
  PolySIM::ThetaL (void)
  /*!
    For simulated data, where 0 is ancenstral, 1 derived.\n
    A simpler version of PolySIM::ThetaL()
    @author Joshua Shapiro
  */
  {
    if(rep->_NumPoly ==0) return 0.;//strtod("NAN",NULL);
    const char state = '1';
    unsigned site;
    unsigned seq;
    int num_changes;
    double thetal = 0.0, nc, nsam = double(rep->_nsam);
    
    for (site = 0; site < rep->_data->numsites (); ++site)
      {
	for (seq = 0, num_changes = 0; seq < unsigned (nsam); ++seq)
	  {
	    num_changes += (*rep->_data)[seq][site] == state ? 1 : 0;
	  }
	nc = double (num_changes);
	thetal +=  nc / (nsam - 1.0);
      }
    return thetal;
  }

  double
  PolySIM::TajimasD (void)
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double Pi = ThetaPi ();
    double W = ThetaW ();
    return ((Pi - W) / Dnominator ());
  }

  double
  PolySIM::Hprime (bool likeThorntonAndolfatto)
  /*!
    Redefinition of PolySNP::Hprime
    @author Joshua Shapiro
  */
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double Hpr = 0.0;
    double pi = ThetaPi ();
    double theta = ThetaW();
    double thetal = ThetaL();
    double a = a_sub_n ();
    double b = b_sub_n ();
    double b_n_plus1 = b_sub_n_plus1();
    double S = rep->_NumPoly;
    double thetasq = (likeThorntonAndolfatto == false) ? S * (S-1)/(a*a + b) : theta*theta;
		  
    double vThetal = 
      (rep->_nsam * theta)/(2.0 * (rep->_nsam - 1.0)) 
      + (2.0 * pow(rep->_nsam/(rep->_nsam - 1.0), 2.0) * (b_n_plus1 - 1.0) - 1.0) * thetasq;
		  
    double vPi = 
      (3.0 * rep->_nsam *(rep->_nsam + 1.0) * theta 
       + 2.0 * ( rep->_nsam * rep->_nsam + rep->_nsam + 3.0) * thetasq  )
      / (9 * rep->_nsam * (rep->_nsam -1.0));
		  
    double cov = 
      (rep->_nsam + 1.0) * theta / (3.0 * (rep->_nsam - 1.0))
      + ( 7.0 * rep->_nsam * rep->_nsam + 3.0 * rep->_nsam - 2.0 - 4.0 * 
	  rep->_nsam *( rep->_nsam + 1.0) * b_n_plus1)
      * thetasq / (2.0 * pow ((rep->_nsam - 1.0), 2.0));
		  
    Hpr = pi - thetal;
    Hpr /= pow ( (vThetal + vPi - 2.0 * cov), 0.5);
    return (Hpr); 
		  
  }

  double
  PolySIM::Dnominator (void)
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double S = rep->_NumPoly;
    double a1, a2, b1, b2, c1, c2, e1, e2;

    a1 = a_sub_n ();
    a2 = b_sub_n ();
    b1 = (rep->_nsam + 1.0) / (3.0 * (rep->_nsam - 1.0));
    b2 = (2.0 * (pow (rep->_nsam, 2.0) + rep->_nsam + 3.0)) / (9.0 * rep->_nsam *
							       (rep->_nsam - 1.0));
    c1 = b1 - 1.0 / a1;
    c2 = b2 - (rep->_nsam + 2.0) / (a1 * rep->_nsam) + a2 / pow (a1, 2.0);
    e1 = c1 / a1;
    e2 = c2 / (pow (a1, 2.0) + a2);
    double denominator = pow ((e1 * S + e2 * S * (S - 1.0)), 0.5);
    return (denominator);
  }


  int
  PolySIM::HudsonsHaplotypeTest (int subsize, int subss)
  /*!
    From Hudson et al (1994) on polymorphism at sod. For simulated data only.
    The function returns a 1 if the number of polymorphisms in a randomly generated subsample
    of the data is less than or equal to subss, 0 otherwise.
    \param subsize the size of the subsample
    \param subss the number of segregating sites in the subsample
    @author Dick Hudson
    @author Kevin Thornton
  */
  {
    int *subslist, i, npoly, seq;
    bool sflag, flag;
    flag = 0;
    sflag = 1;

    subslist = new int[subsize];

    for (i = 0; i < subsize; ++i)
      subslist[i] = i;

    while (sflag)
      {
        npoly = poly (subslist, rep->_nsites,
                      subsize, subss, &seq);
        if (npoly > subss)
          sflag = nextsample (subslist, subsize, int(rep->_nsam), seq);
        else
          {
            flag = 1;
            break;
          }
      }
    if (flag == 1)
      {
        delete[]subslist;
        return (1);
      }
    else
      {
        delete[]subslist;
        return (0);
      }
  }

  int
  PolySIM::poly (int *subslist, int ss,
                 int subsize, int subss, int *seq)
  /*!
    Count the number of polymorphisms in a sample of size subsize.
    Part of the Hudson Haplotype Test. Called by 
    PolySIM::HudsonsHaplotypeTest (const Sequence::SimData * data, int subsize, int subss)
    @author Dick Hudson
    @author Kevin Thornton
  */
  {
    int npoly = 0, i, j, *polyvector;

    polyvector = new int[ss];
    for (i = 0; i < ss; ++i)
      polyvector[i] = 0;

    for (i = 1; i < subsize; ++i)
      for (j = 0; j < ss; ++j)
        if ((*rep->_data)[subslist[i]][j] !=
            (*rep->_data)[subslist[0]][j])
          {
            if (polyvector[j] == 0)
              {
                polyvector[j] = 1;
                npoly++;
                if (npoly == subss + 1)
                  *seq = i;
              }
          }
    delete [] polyvector;
    return (npoly);
  }

  int
  PolySIM::nextsample (int *subslist, int subsize, int nsam, int seq)
  /*!
    Get the next subsample for the HHT.
    Part of the Hudson Haplotype Test.
    Called by PolySIM::HudsonsHaplotypeTest (const Sequence::SimData * data, int subsize, int subss)
    @author Dick Hudson
    @author Kevin Thornton
  */
  {
    int i;
    subslist[seq]++;
    if (subslist[seq] > nsam - subsize + seq)
      {
        seq--;

        if (seq < 0)
          return (0);

        return (nextsample (subslist, subsize, nsam, seq));
      }

    for (i = seq + 1; i < subsize; i++)
      subslist[i] = subslist[i - 1] + 1;
    return (1);
  }

  double
  PolySIM::FuLiD (void)
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double D = 0.0;
    double ExternalMutations = double (NumExternalMutations ());
    double NumMut = double (NumMutations ());
    double a = a_sub_n ();
    double b = b_sub_n ();
    double c = c_sub_n ();

    double vD = 1.0 +
      (pow (a, 2.0) / (b + pow (a, 2.0)) *
       (c - (rep->_nsam + 1.0) / (rep->_nsam - 1.0)));
    double uD = a - 1.0 - vD;
    D = NumMut - a * double (ExternalMutations);
    D /= pow ((uD * NumMut + vD * pow (NumMut, 2.0)), 0.5);
    return (D);
  }

  double
  PolySIM::FuLiF (void)
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double F = 0.0;
    double Pi = ThetaPi ();
    double NumMut = double (NumMutations ());
    double ExternalMutations = double (NumExternalMutations ());
    double a = a_sub_n ();
    double a_n_plus1 = a_sub_n_plus1 ();
    double b = b_sub_n ();
    double c = c_sub_n ();
    double vF = c + 2.0 * (pow (rep->_nsam, 2.0) + rep->_nsam +
                           3.0) / (9.0 * rep->_nsam * (double (rep->_nsam - 1.0)));
    vF -= (2.0 / (rep->_nsam - 1.0));
    vF /= (pow (a, 2.0) + b);

    double uF = 1.0 + (rep->_nsam + 1.0) / (3.0 * (double (rep->_nsam - 1.0)));
    uF -= 4.0 * ((rep->_nsam + 1.0) / (pow (rep->_nsam - 1.0, 2.0))) *
      (a_n_plus1 - 2.0 * rep->_nsam / (rep->_nsam + 1.0));
    uF /= a;
    uF -= vF;

    F = Pi - ExternalMutations;
    F /= pow (uF * NumMut + vF * pow (NumMut, 2.0), 0.5);
    return (F);
  }

  double
  PolySIM::FuLiDStar (void)
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double DStar = 0.0;
    double Singletons = double (NumSingletons ());
    double NumMut = double (NumMutations ());

    double a = a_sub_n ();
    double b = b_sub_n ();
    double d = d_sub_n ();

    double vD = pow (rep->_nsam / (rep->_nsam - 1.0), 2.0) * b;
    vD += pow (a, 2.0) * d;
    vD -= 2.0 * (rep->_nsam * a * (a + 1.0)) /
      (pow (double (rep->_nsam - 1.0), 2.0));
    vD /= (pow (a, 2.0) + b);

    double uD =
      (rep->_nsam / (rep->_nsam - 1.0)) * (a -
                                           (rep->_nsam / (rep->_nsam - 1.0))) - vD;

    DStar = (rep->_nsam / (rep->_nsam - 1.0)) * NumMut - a * double (Singletons);
    DStar /= pow (uD * NumMut + vD * pow (NumMut, 2.0), 0.5);
    return (DStar);
  }

  double
  PolySIM::FuLiFStar (void)
  {
    if(rep->_NumPoly ==0) return strtod("NAN",NULL);
    double FStar = 0.0;
    double Singletons = double (NumSingletons ());
    double Pi = ThetaPi ();
    double NumMut = double (NumMutations ());

    double a = a_sub_n ();
    double a_n_plus1 = a_sub_n_plus1 ();
    double b = b_sub_n ();
    //vF is taken from the correction published by
    //Simonsen et al.  (1995) Genetics 141: 413, eqn A5
    double vF = 2.0 * pow (rep->_nsam, 3.0) + 110.0 * pow (rep->_nsam,  2.0) -  255.0 * rep->_nsam + 153.0;
    vF /= (9.0 * pow (rep->_nsam, 2.0) * (rep->_nsam - 1.0));
    vF += (((2.0 * (rep->_nsam - 1.0) * a) / pow (rep->_nsam, 2.0)) -  (8.0 * b / rep->_nsam));
    vF /= (pow (a, 2.0) + b);

    double uF =  (4.0 * pow (rep->_nsam, 2.0) + 19.0 * rep->_nsam + 3.0 -  12.0 * (rep->_nsam + 1.0) * a_n_plus1);
    uF /= (3.0 * rep->_nsam * (rep->_nsam - 1.0));
    uF /= a;
    uF -= vF;
    FStar = Pi - (((rep->_nsam - 1.0) / rep->_nsam)) * double (Singletons);
    FStar /= pow ((uF * NumMut + vF * pow (NumMut, 2.0)), 0.5);
    return (FStar);
  }

  unsigned
  PolySIM::NumMutations (void)
  /*!
    \return number of mutations in the sample
  */
  {
    return rep->_NumPoly;
  }

  //count the number of singletons
  //only works for infinite sites
  unsigned
  PolySIM::NumSingletons (void)
  /*!
    A version optimized for simulated data where
    character states take on the values 0 or 1.
    \return number of sites where there is a mutation at frequency 1 in the sample
  */
  {
    if (rep->_counted_singletons == false)
      {
        unsigned singletonCount=0;//, changes, j;
	//        unsigned i;
        singletonCount = 0;
	for(std::vector<stateCounter>::const_iterator i = rep->_counts.begin() ;
	    i != rep->_counts.end() ; ++i)
	  {
	    singletonCount += (i->one == 1 || i->zero == 1) ? 1 : 0;
	  }
        rep->_counted_singletons = true;
        rep->_singletons = singletonCount;
        return (singletonCount);
      }
    else
      return rep->_singletons;

    return rep->_singletons;
  }


  unsigned
  PolySIM::NumExternalMutations (void)
  /*!
    similar to num singletons, but it assumes strict
    ancestral vs. derived in the data->
    i.e. for the infinite-sites case, 0 is ancestral,
    1 is derived (as in the case of coalescent simulations).
    note that Sequence does put 1 as the derived state, if 
    you have an outgroup, so this is the routine to use.
    \return the number of derived alleles at frequency 1
  */
  {
    unsigned numExternal = 0;
    for( std::vector<stateCounter>::const_iterator i = rep->_counts.begin() ;
	 i != rep->_counts.end() ; ++i )
      {
	if (i->one == 1)
	  {
	    ++numExternal;
	  }
      }
    return (numExternal);
  }

  unsigned
  PolySIM::Minrec (void)
  /*!
    \return the minimum number of recombination events observed
    in the sample (Hudson and Kaplan 1985). Will return SEQMAXUNSIGNED 
    if there are < 2 segregating sites. 
    \note code is a modification of that provided by Jeff Wall
  */
  {
    if(rep->_NumPoly < 2) return SEQMAXUNSIGNED;
    unsigned a, b, e, numgametes, Rmin=0,x=0;
    bool flag = false;
      
    if (rep->_NumPoly < 2 || x >= (rep->_NumPoly - 1))
      return (0);
    for (a = x + 1; a < rep->_nsites; ++a)
      {
	for (b = (flag == false) ? x : a-1 ; b < a; ++b)
	  {
	    flag = false;
	    numgametes = 0;
	    for (e = 0; e < rep->_nsam; ++e)
	      {
		if ((*rep->_data)[e][b] == '0' &&(*rep->_data)[e][a] == '0')
		  {
		    ++numgametes;
		    break;
		  }
	      }
	    for (e = 0; e < rep->_nsam; ++e)
	      {
		if ((*rep->_data)[e][b] == '0' &&(*rep->_data)[e][a] == '1')
		  {
		    ++numgametes;
		    break;
		  }
	      }
	    for (e = 0; e < rep->_nsam; ++e)
	      {
		if ((*rep->_data)[e][b] == '1' &&(*rep->_data)[e][a] == '0')
		  {
		    ++numgametes;
		    break;
		  }
	      }
	    for (e = 0; e < rep->_nsam; ++e)
	      {
		if ((*rep->_data)[e][b] == '1' &&(*rep->_data)[e][a] == '1')
		  {
		    ++numgametes;
		    break;
		  }
	      }
	    if (numgametes == 4)
	      {
		++Rmin;
		flag = true;
		break;
	      }      
	  }
	if (flag == true)
	  {
	    x = a;
	  }
      }
    return Rmin;
  }

  void PolySIM::WallStats(void)
  {
    unsigned n00 ,n01 ,n10, n11;
    unsigned nhap_curr, nhap_left;
    unsigned n0site1,n0site2;
    nhap_left = SEQMAXUNSIGNED;

    unsigned A = 0;//number of partitions with D' = 1 (see Wall 1999)
    unsigned S = rep->_NumPoly;
    if (S > 1)
      {
	for (unsigned site1 = 0 ; site1 < rep->_nsites-1 ; ++site1)
	  //iterate over sites (actually, adjacent pairs of sites)
	  {
	    for (unsigned site2 = site1+1 ; site2 < rep->_nsites ; ++site2)
	      {
		nhap_curr =0;
		n00 = n01 = n10 = n11 = n0site1 = n0site2 = 0;
		for (unsigned i = 0 ; i < rep->_nsam ; ++i)
		  {
		    switch ( (*rep->_data)[i][site1] )
		      {
		      case '0':
			++n0site1;
			switch ( (*rep->_data)[i][site2] )
			  {
			  case '0':
			    ++n0site2;
			    ++n00;
			    break;
			  case '1':
			    ++n01;
			    break;
			  }
			break;
		      case '1':
			switch ( (*rep->_data)[i][site2] )
			  {
			  case '0':
			    ++n0site2;
			    ++n10;
			    break;
			  case '1':
			    ++n11;
			    break;
			  }
			break;
		      }
		  }
		//the if statement checks to make
		//sure that both sites are polymorphic
		if ( (n0site1 > 0 && n0site1 < rep->_nsam) &&
		     (n0site2 > 0 && n0site2 < rep->_nsam) )
		  {
		    nhap_curr += (n00 > 0) ? 1 : 0;
		    nhap_curr += (n01 > 0) ? 1 : 0;
		    nhap_curr += (n10 > 0) ? 1 : 0;
		    nhap_curr += (n11 > 0) ? 1 : 0;
		    if (site1 == 0)
		      {
			if (nhap_curr == 2)
			  {
			    ++rep->_walls_Bprime;
			    ++A;
			  }
		      }
		    else
		      {
			if (nhap_curr == 2)
			  ++rep->_walls_Bprime;
			if (nhap_curr == 2 && nhap_left != 2)
			  ++A;
		      }
		    site1=site2;
		  }
		nhap_left = nhap_curr;
	      }
	  }
	rep->_walls_B = double(rep->_walls_Bprime)/(double(S-1));
	rep->_walls_Q = (rep->_walls_Bprime + double(A))/(double(S));
      }
    else
      {
	rep->_walls_B = strtod("NAN",NULL);
	rep->_walls_Bprime=0;
	rep->_walls_Q = strtod("NAN",NULL);
      }
    rep->_calculated_wall_stats=true;
  }

  double PolySIM::WallsB(void)
  {
    if (rep->_calculated_wall_stats == false)
      WallStats();
    return rep->_walls_B;
  }
  unsigned PolySIM::WallsBprime(void)
  {
    if (rep->_calculated_wall_stats == false)
      WallStats();
    return rep->_walls_Bprime;
  }
  double PolySIM::WallsQ(void)
  {
    if (rep->_calculated_wall_stats == false)
      WallStats();
    return rep->_walls_Q;
  }
}
