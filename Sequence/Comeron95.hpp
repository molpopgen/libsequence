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

#ifndef COMERON95_H
#define COMERON95_H
/*! \file Comeron95.hpp
  @brief declaration of Sequence::Comeron95 to calculate Ka/Ks
*/
/*!
  \class Sequence::Comeron95 Sequence/Comeron95.hpp
  \ingroup kaks
  A definition of an object to implement Comeron's (1995) method to calculate Ka and Ks
  The reference is:
  Comeron, J. (1995) J. Molecular Evolution, 41: 1152-1159.  The object has a simple purpose,
  but it requires a lot of code to do it.  The calculation depends upon the following classes:\n
  Sequence::RedundancyCom95\n
  Sequence::Sites\n
  Sequence::SingleSub\n
  Sequence::TwoSubs\n
  Sequence::ThreeSubs\n
  \n
  The calculation can be done on any objects of type Sequence::Seq *, or anything derived
  from there.\n
  \n
  You should be warned that this class depends heavily on namespace Sequence (the libsequence
  library), and may throw any and all exceptions defined in that library, so you should see the
  documentation for that library.
  \note This class allows weighting of pathways where codons differ at more than 1 site to be
  done using any arbitrary weighting scheme.  To do so, however, requires deriving classes
  from Sequence::WeightingScheme2 and Sequence::WeightingScheme2.  See the code for 
  Sequence::GranthamWeights2 and Sequence::GranthamWeights3 for an example of how to do this 
  (in the files GranthamWeights.cc and GranthamWeights.hpp).  You will also need to study
  TwoSubs.cc/.hpp and ThreeSubs.cc/.hpp to understand how branches and intermediate codons
  are labelled (because you need to do it correctly in your implementation of weighting 
  schemes)\n
  \note This class is best used in the form \c auto_ptr<Comeron95>, as illustrated
  in the example below.  The reason to use \c auto_ptr is to prevent resource leaks in
  case an exception is thrown.  This class has strict requirements about the data it
  is passed, resulting in several places where exceptions can be thrown.  Most of
  these situations are caught in the constructor, which will prevent most leaks.
  However, anything is possible, so use \c auto_ptr
  \n
  Example: a simple program to calculate Ka and Ks:\n
  \code
  #include<iostream>
  #include<memory>
  #include <Sequence/Fasta.hpp>
  #include <Sequence/Alignment.hpp>
  #include <Sequence/SeqExceptions.hpp>
  #include <Sequence/Comeron95.hpp>
 
  //using namespace std;
  using namespace Sequence;
 
  int main
  {
    char *filename = argv[1];
    vector<Fasta > data;
    
    //read data in from a file in Fasta Format
    //program will exit if file is badly format
    //the reason for the exit will be an uncaught
    //exception Sequence::badFormat
    Alignment::GetData(data,filename);
 
    //abort() unless all sequence objects are the same length
    assert(Alignment::IsAlignment(data));
 
    //iterate over the data
    for(unsigned i = 0 ; i < data.size()-1 ; ++i)
    {
      for(unsigned j = i+1 ; j < data.size() ; ++j)
      {
      //calculate and output results
	try
	{
        auto_ptr<Comeron95> C(new Comeron95(&data[i],&data[j]));
	cout << data[i].GetName() << '\t' << data[j].GetName() << '\t';
	cout << C->ka() << '\t' << C->ks() << '\t' << C->ratio() << endl;
	if (C->ratio() > 1.0 && C->ratio() != 999.0)
	  cout << "congratulations, you win!" << endl;
	}
	catch (SeqException &e)
	{
	  cout << e << endl;
	}
      }
    }
  }
  \endcode
  @short Ka and Ks by Comeron's (1995) method
*/

/*! \example gestimator.cc
 */
#include <Sequence/SeqEnums.hpp>
#include <boost/utility.hpp>

namespace Sequence
  {
  class Seq;
  class Sites;
  class RedundancyCom95;
  class WeightingScheme2;
  class WeightingScheme3;
    class Comeron95 : boost::noncopyable
    {
    private:
      bool __2wasNULL,__3wasNULL,__red_was_NULL;
      WeightingScheme2 *weights2;
      WeightingScheme3 *weights3;
      int valid, maxhits, genetic_code, weighting_scheme;
      //see Comeron '95 for a discussion of the method--will document later
      double Qs, Bs, Qa, Ba, A2S, A4, As, A2V, A0, Aa;
      double q0, q2S, q2V, q4, p0, p2S, p2V, p4;
      Sites *sites;
      const RedundancyCom95 *sitesObj;
      void diverge (const Sequence::Seq * seq1, const Sequence::Seq * seq2,
                    WeightingScheme2 *_weights2,
                    WeightingScheme3 *_weights3);
      void omega (const Sequence::Seq * seqobj1, const Sequence::Seq * seqobj2);
      double Ka, Ks;
    public:
      explicit Comeron95 (const Sequence::Seq * seqa,
                          const Sequence::Seq * seqb,
                          int max = 3, 
			  const Sequence::RedundancyCom95 * genetic_code_redundancy = NULL,
			  GeneticCodes code = UNIVERSAL,
                          WeightingScheme2 *weights2 = NULL,
                          WeightingScheme3 *weights3 = NULL);
      ~Comeron95 (void);
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
      double L0 (void) const;
      double L2S (void) const;
      double L2V (void) const;
      double L4 (void) const;
  };
}

#endif
