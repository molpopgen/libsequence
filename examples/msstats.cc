/* 
   msstats - read data from ms via stdin, calculate common summary statistics

   Copyright (C) 2002 Kevin Thornton

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

*/

#include <iostream>
#include <vector>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <cstdio>

using namespace std;
using namespace Sequence;

int main(int argc, char *argv[]) 
{
  SimParams p;
  cin >> p;
  SimData d(p.totsam());

  std::ios_base::sync_with_stdio(true);

  int rv;
  while( (rv = d.fromfile(stdin)) != EOF )
    {
      PolySIM P(&d);
      cout <<P.NumPoly()  << '\t' 
	   << P.ThetaW()  << '\t' 
	   << P.ThetaPi() << '\t'
	   << P.ThetaH()  << '\t' 
	   << P.TajimasD() << '\t'
	   << P.FuLiD() << '\t'
	   << P.FuLiF() << '\t'
	   << P.FuLiDStar() << '\t'
	   << P.FuLiFStar() << endl;
    }
}
