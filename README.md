#libsequence - A C++ class library for evolutionary genetic analysis



  Copyright (C) 2002 Kevin Thornton

  libsequence2 is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Comments are welcome.

	- Kevin Thornton <krthornt@uci.edu>

#Citation

If you use the library for your research, please cite:

@article{libsequence,
author = {Thornton, Kevin},
title = {{Libsequence: a C++ class library for evolutionary genetic analysis.}},
journal = {Bioinformatics (Oxford, England)},
year = {2003},
volume = {19},
number = {17},
pages = {2325--2327},
month = nov
}

The manuscript is available online at http://bioinformatics.oxfordjournals.org/content/19/17/2325.short

#Installation

##Dependencies
GNU Scientific Library: http://gnu.org/software/gsl

Boost C++ libraries: http://www.boost.org

zlib: http://zlib.net

##Simplest installation instructions

The "CXXFLAGS= " part of the command is to over-ride ./configure's desire to add -g -O2 to the compile options

CXXFLAGS= ./configure

CXXFLAGS= make

sudo make install
