# libsequence - A C++ class library for evolutionary genetic analysis



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

# User's group

Please post to the [libsequence user group](https://groups.google.com/forum/#!forum/libsequence-users) for help.

# Build status

* master branch: [![Build Status](https://travis-ci.org/molpopgen/libsequence.svg?branch=master)](https://travis-ci.org/molpopgen/libsequence)
* dev branch: [![Build Status](https://travis-ci.org/molpopgen/libsequence.svg?branch=dev)](https://travis-ci.org/molpopgen/libsequence)


# Citation

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

## Revision history.

The revision history of the library is [here](REVISION_HISTORY.md).  The document describes what changed for a given release.

## Obtaining the source code

### Obtaining the master branch
You have a few options:

* Clone the repo (best option): git clone https://github.com/molpopgen/libsequence.git
* Click on "Download Zip" at https://github.com/molpopgen/libsequence


### Obtaining a specific release
Again, a few options:

* Click on "Releases" at https://github.com/molpopgen/libsequence, then download the one you want
* Clone the repo (see previous section)
* Get a list of releases by saying "git tag -l"
* Checkout the release you want.  For example "git checkout 1.8.0"

## Installation

### Dependencies

1. A C++11-compliant compiler (see next section)
2. zlib: http://zlib.net
3. [Intel's TBB library](http://threadbuildingblocks.org) is used for parallelization.

### Compilers

I support the following compilers:

* [GCC](http://gcc.gnu.org)
* [clang](http://clang.llvm.org)

The following compilers are not supported:

* [Intel](https://software.intel.com/en-us/intel-compilers)

The Intel compiler suffers from the following issues:

* It appears to no longer be free for academic use. (boo!)
* It appears to be based on a version of libstdc++ that is too old to be compatible with libsequence (see [here](https://github.com/molpopgen/libsequence/pull/4) for some discussion of the issue)

### Optional dependencies

1. [htslib](http://htslib.org) The configure script will attempt to detect the presence of htslib on your system.  If the library is present, then libsequence will compile with support for features like direct reading from BAM files.  If htslib is not present, those features will not be compiled.
   
### Simplest installation instructions

~~~
./configure
make
sudo make install
~~~

The build conditions can be adjusted via the usual environment variables.  To compile an optimized "release" build:

~~~
./configure CXXFLAGS="-O3 -DNDEBUG"
~~~

To compile a debugger-friendly build:

~~~
./configure CXXFLAGS="-O0 -g"
~~~

To change the compiler, set the C and C++ compiler variables:

~~~
./configure CC=gcc CXX=g++
~~~

#### Compiling unit tests and examples

To compile unit testing suite and example programs

~~~
make check
~~~

or

~~~
cd test
make check
~~~

Note that the library must be built prior to "make check", but _you do not have to install the library prior ot "make check"_.  The examples and unit tests are statically-linked to the version of the library that will be found in src/.libs after a "make" command.  I do this so that one can perform unit tests without having to install the library.  I use static linking here to avoid any possible confusion with an existing libsequence installation.

#### Running the unit tests

~~~
cd test && sh runTests.sh
~~~

### More complex installation scenarios

Some users may not have the dependent libraries installed in the standard locations on their systems.  Note that "standard" means wherever the compiler system looks for header files during compilation and libraries during linking.  This scenario is common on OS X systems where users have used some sort of "system" to install various libraries rather than installing from source directly.  In order to accomodate such situations, the user must provide the correct path to the include and lib directories.  For example, assume that the dependend libraries are in /opt on your system.  You would install libsequence as follows:

CPPFLAGS=-I/opt/include LDFLAGS="$LDFLAGS -l/opt/lib" ./configure

make

~~~
sudo make install
~~~

Note that the modification of LDFLAGS prepends the current value of LDFLAGS if it exists.  This allows for scenarios where the system's search path for libraries may have been modified by the user or sysadmin via a modification of that shell variable.  (One could also do the same with CPPFLAGS, FYI.)

### Installing libsequence locally

If you do not have permission to "sudo make install", you can install the library in your $HOME:

./configure --prefix=$HOME

Then, when compiling any program using libsequence, you need to add

~~~
-I$HOME/include
~~~
to any compilation commands and

~~~
-L$HOME/lib
~~~

to any linking commands.

When running programs linking to any of the above run-time libraries, and depending on your system, you may also need to adjust variables like LD_LIBRARY_PATH to prepend $HOME/lib to them, etc., but you'll need to figure that out on case-by-case basis, as different systems can behave quite differently.

### Another installation option (not supported by the libsequence author)

I've recently been made aware that there is a method for installing libsequence using the [brew.sh](http://brew.sh/) system.  This system allows the [homebrew-science](https://github.com/Homebrew/homebrew-science) git repo to be used to obtain libsequence.  I do not use this system myself, nor do I know how to.

## Using libsequence to compile other programs

If libsequence is not installed in a standard path, then you must provide the appropriate include (-I) and link path (-L) commands to your compiler.  This may be done in various ways, e.g., via a configure script or your own Makefile.

A program that depends on libsequence must provide at least the following libraries to the linker:

~~~
-lsequence -lz 
~~~

If you are using features depending on htslib, the linking options become

~~~
-lsequence -lz -lhts
~~~

