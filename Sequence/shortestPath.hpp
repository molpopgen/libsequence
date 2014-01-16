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

#ifndef __SHORTEST_PATH_HPP__
#define __SHORTEST_PATH_HPP__
#include <vector>
#include <utility>
#include <memory>
#include <Sequence/SeqEnums.hpp>
#include <Sequence/SeqExceptions.hpp>
#include <boost/tuple/tuple.hpp>
/*! \file shortestPath.hpp
  @short Routines to find the shortest distance between any 2 codons, using
  Grantham's distance. Declares the class Sequence::shortestPath, and the
  functions Sequence::mutsShortestPath and Sequence::diffType
*/

/*!
  \class Sequence::shortestPath Sequence/shortestPath.hpp
  A class which calculates the shortest path between two
  codons.  The length of a path is in terms of the sum
  of the Grantham's distances along it's branches.
  \note It is often the case (esp. when codons differ at
  all three sites) that the various paths are of equal
  length.  What this means is that the order in which
  the changes occur don't matter.  In such cases,
  one of the paths is returned arbitrarily.
  \short Calculate shortest path between 2 codons
  \ingroup CodonPaths
*/
namespace Sequence
{
  class shortestPathImpl;
  class shortestPath 
  {
  private:
    std::auto_ptr<shortestPathImpl> impl;
  public:
    /*!
      An enum type to describe the shortest path between 2 codons.
      An S refers to a synonymous change, N nonsynonymous. For
      example, SSN means that the shortest path b/w 2 codons
      requires 2 synonymous and 1 nonsynonymous change. If the
      2 codons don't differ, the value is NONE.  If the path
      cannot be determined, the value is AMBIG
    */
    enum pathType {S,N,SS,SN,NN,SSS,SSN,SNN,NNN,NONE,AMBIG};

    explicit shortestPath(const std::string &codon1,
			  const std::string &codon2,
			  const Sequence::GeneticCodes & code = Sequence::UNIVERSAL) 
      ;
    ~shortestPath();
    pathType type() const;
    double path_distance() const;
    typedef std::vector<std::string>::const_iterator const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
  };

  std::pair<unsigned,unsigned> mutsShortestPath(const std::string &codon1,
						const std::string &codon2,
						const Sequence::GeneticCodes 
						& code = Sequence::UNIVERSAL)
    ;

  std::pair<unsigned,shortestPath::pathType> diffType(const std::string &codon1,
						      const std::string &codon2,
						      const Sequence::GeneticCodes 
						      & code = Sequence::UNIVERSAL)
    ;

  boost::tuple<shortestPath::pathType,shortestPath::pathType,shortestPath::pathType>
  diffTypeMulti(const std::string &codon1,
		const std::string &codon2,
		const Sequence::GeneticCodes 
		& code = Sequence::UNIVERSAL)
    ;
}
#endif
