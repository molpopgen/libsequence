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
#include <Sequence/Grantham.hpp>
#include <Sequence/Translate.hpp>
#include <Sequence/PathwayHelper.hpp>
#include <Sequence/GranthamWeights.hpp>

using std::string;

namespace Sequence
  {
  GranthamWeights2::GranthamWeights2(GeneticCodes genetic_code)
      :code(genetic_code)
      /*!
        \param genetic_code the code used to Translate the codons
      */
  {
  }

  GranthamWeights2::~GranthamWeights2(void)
  {}

  void GranthamWeights2::Calculate(const std::string &codon1, const std::string &codon2) const
  /*!
    Calculate actually calculates the weights for each branch
    \param codon1 a std::string of length 3 representing a sense codon
    \param codon2 a std::string of length 3 representing a sense codon
  */
    {
      Grantham gdist;
      string intermediates[2];
      Intermediates2(intermediates,codon1,codon2);
      //assign weights to pathways
      //weights are assiged by the total length of each path, as
      //measured by Grantham distances
      double len_path_1 = 0.0, len_path_2 = 0.0;

      string t1 = Translate (codon1.begin(),codon1.end(), code);
      string t2 = Translate (intermediates[0].begin(),
                             intermediates[0].end(), code);
      len_path_1 += gdist ((t1)[0], (t2)[0]);

      t1 = Translate (intermediates[0].begin(),intermediates[0].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      len_path_1 += gdist ((t1)[0], (t2)[0]);

      t1 = Translate (codon1.begin(),codon1.end(), code);
      t2 = Translate (intermediates[1].begin(),intermediates[1].end(), code);
      len_path_2 += gdist ((t1)[0], (t2)[0]);

      t1 = Translate (intermediates[1].begin(),intermediates[1].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      len_path_2 += gdist ((t1)[0], (t2)[0]);

      //calculate the weights themselves
      //double w_path1 = 0., w_path2 = 0.,
      double w_tot = 0.;
      __weights[0]=0.;
      __weights[1]=0.;

      if (fabs(len_path_1-0.) <= DBL_EPSILON && fabs(len_path_2-0.) <= DBL_EPSILON)
        {
          __weights[0] = __weights[1] = 0.5;
          w_tot = 1.;
        }
      else
        {
          __weights[0] = 1. - len_path_1 / (len_path_1 + len_path_2);
          __weights[1] = 1. - len_path_2 / (len_path_1 + len_path_2);
          w_tot = __weights[0] + __weights[1];
        }

      __weights[0] /= w_tot;
      __weights[1] /= w_tot;

    }

  double * GranthamWeights2::weights(void) const
  /*!
    \return a double * of size 2 (1 value for each branch)
  */
    {
      return __weights;
    }

  GranthamWeights3::GranthamWeights3(GeneticCodes genetic_code)
      :code(genetic_code)
      /*!
        \param genetic_code the code used to Translate the codons
      */
  {
  }

  GranthamWeights3::~GranthamWeights3(void)
  {}

  void GranthamWeights3::Calculate(const std::string &codon1, const std::string &codon2) const
  /*!
    Calculate actually calculates the weights for each branch
    \param codon1 a std::string of length 3 representing a sense codon
    \param codon2 a std::string of length 3 representing a sense codon
  */
    {
      Grantham gdist;
      string intermediates[9];
      Intermediates3(intermediates,codon1,codon2);
      double len_path_1 = 0.0, len_path_2 = 0.0, len_path_3 =
                                              0., len_path_4 = 0., len_path_5 = 0., len_path_6 = 0.;
      double dist = 0.;

      string t1 = Translate (codon1.begin(), codon1.end(),code);
      string t2 = Translate (intermediates[0].begin(), intermediates[0].end(),code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_1 += dist;


      t1 = Translate (intermediates[0].begin(), intermediates[0].end(), code);
      t2 = Translate (intermediates[1].begin(), intermediates[1].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_1 += dist;


      t1 = Translate (intermediates[1].begin(), intermediates[1].end(),code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_1 += dist;

      //path2

      t1 = Translate (codon1.begin(),codon1.end(), code);
      t2 = Translate (intermediates[0].begin(),intermediates[0].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_2 += dist;

      t1 = Translate (intermediates[0].begin(),intermediates[0].end(), code);
      t2 = Translate (intermediates[2].begin(),intermediates[2].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_2 += dist;

      t1 = Translate (intermediates[2].begin(),intermediates[2].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_2 += dist;

      //path 3
      t1 = Translate (codon1.begin(),codon1.end(), code);
      t2 = Translate (intermediates[3].begin(),intermediates[3].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_3 += dist;

      t1 = Translate (intermediates[3].begin(), intermediates[3].end(),code);
      t2 = Translate (intermediates[4].begin(), intermediates[4].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_3 += dist;

      t1 = Translate (intermediates[4].begin(), intermediates[4].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_3 += dist;

      //path4
      t1 = Translate (codon1.begin(),codon1.end(), code);
      t2 = Translate (intermediates[3].begin(), intermediates[3].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_4 += dist;

      t1 = Translate (intermediates[3].begin(), intermediates[3].end(), code);
      t2 = Translate (intermediates[5].begin(), intermediates[5].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_4 += dist;

      t1 = Translate (intermediates[5].begin(), intermediates[5].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_4 += dist;

      //path 5
      t1 = Translate (codon1.begin(),codon1.end(), code);
      t2 = Translate (intermediates[6].begin(), intermediates[6].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_5 += dist;

      t1 = Translate (intermediates[6].begin(), intermediates[6].end(), code);
      t2 = Translate (intermediates[7].begin(), intermediates[7].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_5 += dist;

      t1 = Translate (intermediates[7].begin(), intermediates[7].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_5 += dist;

      //path 6
      t1 = Translate (codon1.begin(),codon1.end(), code);
      t2 = Translate (intermediates[6].begin(), intermediates[6].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_6 += dist;

      t1 = Translate (intermediates[6].begin(), intermediates[6].end(), code);
      t2 = Translate (intermediates[8].begin(), intermediates[8].end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_6 += dist;

      t1 = Translate (intermediates[8].begin(), intermediates[8].end(), code);
      t2 = Translate (codon2.begin(),codon2.end(), code);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_6 += dist;

      __weights[0] = 1. / len_path_1;
      __weights[1] = 1. /len_path_2;
      __weights[2] = 1. /len_path_3;
      __weights[3] = 1. / len_path_4;
      __weights[4] = 1. /len_path_5;
      __weights[5] = 1./ len_path_6;

      //scale weights to sum to 1
      double w_tot =
        __weights[0] + __weights[1] + __weights[2] + __weights[3] + __weights[4] + __weights[5];
      __weights[0] /= w_tot;
      __weights[1] /= w_tot;
      __weights[2] /= w_tot;
      __weights[3] /= w_tot;
      __weights[4] /= w_tot;
      __weights[5] /= w_tot;
    }

  double * GranthamWeights3::weights(void) const
  /*!
    \return a double * of size 6 (1 value for each branch)
  */
    {
      return __weights;
    }
}
