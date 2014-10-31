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

##Revision history.

The revision history of the library is [here](REVISION_HISTORY.md).  The document describes what changed for a given release.

##Obtaining the source code

###Obtaining the master branch
You have a few options:
<ol>
<li> Clone the repo (best option): git clone https://github.com/molpopgen/libsequence.git</li>
<li> Click on "Download Zip" at https://github.com/molpopgen/libsequence </li>
</ol>

###Obtaining a specific release
Again, a few options:
<ol>
<li> Click on "Releases" at https://github.com/molpopgen/libsequence, then download the one you want </li>
<li> Clone the repo (see previous section)</li>
<ol>
<li> Get a list of releases by saying "git tag -l" </li>
<li> Checkout the release you want.  For example "git checkout 1.8.0"</li>
</ol>
</ol>

##Installation

###Dependencies

1. zlib: http://zlib.net
2. A C++11-compliant compiler

###Optional dependencies

1. [htslib](http://htslib.org) The configure script will attempt to detect the presence of htslib on your system.  If the library is present, then libsequence will compile with support for features like direct reading from BAM files.  If htslib is not present, those features will not be compiled.
   
###Simplest installation instructions

```
./configure
make
sudo make install
```
###More complex installation scenarios

Some users may not have the dependent libraries installed in the standard locations on their systems.  Note that "standard" means wherever the compiler system looks for header files during compilation and libraries during linking.  This scenario is common on OS X systems where users have used some sort of "system" to install various libraries rather than installing from source directly.  In order to accomodate such situations, the user must provide the correct path to the include and lib directories.  For example, assume that the dependend libraries are in /opt on your system.  You would install libsequence as follows:

CXXFLAGS=-I/opt/include LDFLAGS="$LDFLAGS -l/opt/lib" ./configure

make

sudo make install

Note that the modification of LDFLAGS prepends the current value of LDFLAGS if it exists.  This allows for scenarios where the system's search path for libraries may have been modified by the user or sysadmin via a modification of that shell variable.  (One could also do the same with CXXFLAGS, FYI.)

###Installing libsequence locally

If you do not have permission to "sudo make install", you can install the library in your $HOME:

./configure --prefix=$HOME

Then, when compiling any program using libsequence, gsl, and/or boost, you need to add

> -I$HOME/include

to any compilation commands and

> -L$HOME/lib

to any linking commands.

When running programs linking to any of the above run-time libraries, and depending on your system, you may also need to adjust variables like LD _ LIBRARY _ PATH to prepend $HOME/lib to them, etc., but you'll need to figure that out on case-by-case basis, as different systems can behave quite differently.

###Another installation option (not supported by the libsequence author)

I've recently been made aware that there is a method for installing libsequence using the [brew.sh](http://brew.sh/) system.  This system allows the [homebrew-science](https://github.com/Homebrew/homebrew-science) git repo to be used to obtain libsequence.  I do not use this system myself, nor do I know how to.

##Using libsequence to compile other programs

If libsequence is not installed in a standard path, then you must provide the appropriate include (-I) and link path (-L) commands to your compiler.  This may be done in various ways, e.g., via a configure script or your own Makefile.

A program that depends on libsequence must provide at least the following libraries to the linker:

-lsequence -lz -lgsl -lgslcblas

#Compiling the examples

There are several example programs in the examples subdirectory.  If you have installed libsequence in a standard path (e.g., /usr/local/lib and /usr/local/include), then you compile the examples by saying

```
make
```

If you have installed the library elsewhere, such as $HOME, then you need to adjust LDFLAGS as follows:

```
LDFLAGS=-L$HOME/lib make
```

If you have dependencies like boost, gsl, in locations other than /usr/local (or their moral equivalent on your system, then you will likely need to manually edit CXXFLAGS in the Makefile to add a -I flag to the folder containing header files.  For example, you may change the variable from

> CXXFLAGS = -O3 -Wall -W -I..

to

> CXXFLAGS = -O3 -Wall -W -I.. -I$(HOME)/include
