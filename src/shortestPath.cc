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

#include <Sequence/shortestPath.hpp>
#include <Sequence/PathwayHelper.hpp>
#include <Sequence/Translate.hpp>
#include <Sequence/Grantham.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/SeqProperties.hpp>
#include <algorithm>

namespace Sequence
{
#ifndef DOXYGEN_SKIP //doxygen should skip this
  struct shortestPathImpl
  {
    shortestPath::pathType _type;
    Sequence::Grantham distances;
    double _distance;
    std::string t1,t2;
    std::vector<std::string> _path;
    std::pair<double,shortestPath::pathType>process_path(const std::string &intermediate,
							 const Sequence::GeneticCodes & code);
    std::pair<double,shortestPath::pathType>process_path2(const std::string &intermediate1,
							  const std::string &intermediate2,
							  const Sequence::GeneticCodes & code);
    shortestPathImpl(const std::string &codon1,
		     const std::string &codon2,
		     const Sequence::GeneticCodes & code);
  };

  shortestPathImpl::shortestPathImpl(const std::string &codon1,
				     const std::string &codon2,
				     const Sequence::GeneticCodes & code):
    _type(shortestPath::AMBIG),
    distances(Sequence::Grantham()),
    _distance(0.)
  {
    if(codon1.length()!=3 || codon2.length() != 3)
      {
	throw Sequence::SeqException("Codons are not both of length 3");
      }
    t1 = Sequence::Translate(codon1.begin(),codon1.end(),code);
    t2 = Sequence::Translate(codon2.begin(),codon2.end(),code);

    _path.push_back(codon1);

    unsigned ndiffs = Sequence::NumDiffs(codon1,codon2);
    if (ndiffs == 0)
      {
	_type = shortestPath::NONE;
      }
    else if (t1 != "X" && t2 != "X" && 
	     (std::find_if(codon1.begin(),codon1.end(),
			   ambiguousNucleotide()) == codon1.end()) &&
	     (std::find_if(codon2.begin(),codon2.end(),
			   ambiguousNucleotide()) == codon2.end()) )
      {
	if (ndiffs == 1)
	  {
	    _distance = distances(t1[0],t2[0]);
	    if(t1 != t2)
	      {
		_type = shortestPath::N;
	      }
	    else
	      {
		_type = shortestPath::S;
	      }

	  }
	else if (ndiffs == 2)
	  {
	    std::string intermediates[2];
	    Sequence::Intermediates2(intermediates,codon1,codon2);
	    std::string t_int1 = Sequence::Translate(intermediates[0].begin(),
						     intermediates[0].end());
	    std::string t_int2 = Sequence::Translate(intermediates[1].begin(),
						     intermediates[1].end());

	    std::vector< std::pair<double, shortestPath::pathType> > 
	      paths(2,std::pair<double, shortestPath::pathType>(0.,shortestPath::AMBIG));

	    //update pathway lengths
	    //process pathway 1 (codon1->intermediates[0]->codon2);
	    paths[0] = process_path(intermediates[0],code);
	    //process pathway 2 (codon1->intermediates[1]->codon2)
	    paths[1] = process_path(intermediates[1],code);
	    if (paths[0].first < paths[1].first)
	      {
		_type = paths[0].second;
		_distance = paths[0].first;
		_path.push_back(intermediates[0]);
	      }
	    else
	      {
		_type = paths[1].second;
		_distance = paths[1].first;
		_path.push_back(intermediates[1]);
	      }
	  }
	else if (ndiffs == 3)
	  {
	    std::string intermediates[9];
	    Sequence::Intermediates3(intermediates,codon1,codon2);
	    std::vector< std::pair<double, shortestPath::pathType> > 
	      paths(6,std::pair<double, shortestPath::pathType>(0.,shortestPath::AMBIG));
	    //there are 6 paths b/w codons that differ at all 3 sites
	    //the indexing of array intermediates is describing in the documentation 
	    //for class Sequence::ThreeSubs in the manual
	    paths[0] = process_path2(intermediates[0],intermediates[1],code);
	    paths[1] = process_path2(intermediates[0],intermediates[2],code);
	    paths[2] = process_path2(intermediates[3],intermediates[4],code);
	    paths[3] = process_path2(intermediates[3],intermediates[5],code);
	    paths[4] = process_path2(intermediates[6],intermediates[7],code);
	    paths[5] = process_path2(intermediates[6],intermediates[8],code);
	      
	    size_t shortest = 0;
	    double shortest_path = paths[0].first;
	    for(unsigned i = 0 ; i < 6 ; ++i)
	      {
		if (paths[i].first < shortest_path)
		  {
		    shortest=i;
		    shortest_path=paths[i].first;
		  }
	      }
	    _type = paths[shortest].second;
	    _distance = paths[shortest].first;
	    switch (shortest)
	      {
	      case 0:
		_path.push_back(intermediates[0]);
		_path.push_back(intermediates[1]);
		break;
	      case 1:
		_path.push_back(intermediates[0]);
		_path.push_back(intermediates[2]);
		break;
	      case 2:
		_path.push_back(intermediates[3]);
		_path.push_back(intermediates[4]);
		break;
	      case 3:
		_path.push_back(intermediates[3]);
		_path.push_back(intermediates[5]);
		break;
	      case 4:
		_path.push_back(intermediates[6]);
		_path.push_back(intermediates[7]);
		break;
	      case 5:
		_path.push_back(intermediates[6]);
		_path.push_back(intermediates[8]);
		break;
	      }
	  }
      }
    else
      {
	_type = shortestPath::AMBIG;
      }
    _path.push_back(codon2);
  }

  std::pair<double,shortestPath::pathType> 
  shortestPathImpl::process_path(const std::string &intermediate,
				 const Sequence::GeneticCodes & code)
  {
    std::string tint = Sequence::Translate(intermediate.begin(),intermediate.end(),
					   code);
    shortestPath::pathType t = shortestPath::AMBIG;
    double d = ( distances(t1[0],tint[0]) + distances(tint[0],t2[0]) );
    if ( (t1 != tint) && (tint != t2) )
      {
	t = shortestPath::NN;
      }
    else if ( ((t1 != tint) && (tint == t2)) ||
	      ((t1 == tint) && (tint != t2)) )
      {
	t = shortestPath::SN;
      }
    else if ( t1==tint && tint==t2 )
      {
	t = shortestPath::SS;
      }
    return std::make_pair(d,t);
  }

  std::pair<double,shortestPath::pathType> 
  shortestPathImpl::process_path2(const std::string &intermediate1,
				  const std::string &intermediate2,
				  const Sequence::GeneticCodes & code)
  {
    std::string tint1 = Sequence::Translate(intermediate1.begin(),
					    intermediate1.end(),code);
    std::string tint2 = Sequence::Translate(intermediate2.begin(),
					    intermediate2.end(),code);
    shortestPath::pathType t = shortestPath::AMBIG;
    double d = (distances(t1[0],tint1[0]) + 
		distances(tint1[0],tint2[0]) + 
		distances(tint2[0],t2[0]));

    bool a = (t1==tint1),
      b = (tint1==tint2),
      c = (tint2==t2);

    if ( !a && !b && !c )
      {
	t = shortestPath::NNN;
      }
    else if ( (a && b && !c) ||
	      (a && !b && c) ||
	      (!a && b && c) )
      {
	t = shortestPath::SSN;
      }
    else if ( (!a && !b && c) ||
	      (!a && b && !c) ||
	      ( a && !b && !c) )
      {
	t = shortestPath::SNN;
      }

    else if (a && b && c)
      {
	t = shortestPath::SSS;
      }
    return std::make_pair(d,t);
  }
#endif

  shortestPath::shortestPath(const std::string &codon1,
			     const std::string &codon2,
			     const Sequence::GeneticCodes & code) 
    

  /*!
    \param codon1 a std::string of length 3
    \param codon2 a std::string of length 3
    \param code which genetic code to use
    \pre (codon1.length() == 3 && codon2.length() ==3)
    \note If either codon1 or codon2 contain characters other than {A,G,C,T},
    the pathway type will be assigned shortestPath::AMBIG
  */
  {
    try
      {
	impl = std::auto_ptr<shortestPathImpl>(new shortestPathImpl(codon1,codon2,code));
      }
    catch(Sequence::SeqException &e)
      {
	throw;
      }
    catch(std::exception &e)
      {
	throw (Sequence::SeqException(e.what()));
      }
    catch (...)
      {
	throw (Sequence::SeqException("caught exception of unknown type"));
      }
  }

  shortestPath::~shortestPath()
  {
  }

  shortestPath::pathType shortestPath::type() const
  /*!
    \returns a value from the enum type shortestPath::pathType
    representing the type of the shortest path.
  */
  { 
    return impl->_type;
  }

  shortestPath::const_iterator shortestPath::begin() const
  /*!
    \return an iterator to the beginning of the shortest path
    between the 2 codons. The value type if the iterator
    is std::string
  */
  {
    return impl->_path.begin();
  }

  shortestPath::const_iterator shortestPath::end() const
  /*!
    \return an iterator to the end of the shortest path
    between the 2 codons. The value type if the iterator
    is std::string
  */
  {
    return impl->_path.end();
  }

  double shortestPath::path_distance() const
  /*!
    \return the total Grantham's distance of the shortest path
  */
  {
    return impl->_distance;
  }

  std::pair<unsigned,unsigned> mutsShortestPath(const std::string &codon1,
						const std::string &codon2, 
						const Sequence::GeneticCodes & code)
    

  /*!
    \return a std::pair<unsigned,unsigned> representing the number of silent
    and replacement changes b/w 2 codons, as calculated by Sequence::shortestPath.
    The first member of the pair is the number of silent changes in the shortest path,
    and the second the number of replacement changes. If the pathway type is
    shortestPath::AMBIG, both members of the return value will be equal to
    Sequence::SEQMAXUNSIGNED, which is declared in <Sequence/SeqConstants.hpp>

    For example:
    \code
    #include <Sequence/shortestPath.hpp>
    #include <iostream>
    int main(int argc, char **argv)
    {
    //the shortest path between AAA and GGG, using the universal code,
    //requires 1 silent and 2 replacement changes.
    std::pair<unsigned,unsigned> muts = Sequence::mutsShortestPath("AAA","GGG");
    std::cout << muts.first
    << ' '
    << muts.second
    << std::endl;
    }
    \endcode
    \ingroup CodonPaths
  */
  {
    try
      {
	shortestPath sp(codon1,codon2,code);
	shortestPath::pathType type = sp.type();
	switch (type)
	  {
	  case shortestPath::S :
	    {
	      return std::make_pair(1,0);
	      break;
	    }
	  case shortestPath::N :
	    {
	      return std::make_pair(0,1);
	      break;
	    }
	  case shortestPath::SS :
	    {
	      return std::make_pair(2,0);
	      break;
	    }
	  case shortestPath::SN :
	    {
	      return std::make_pair(1,1);
	      break;
	    }
	  case shortestPath::NN :
	    {
	      return std::make_pair(0,2);
	      break;
	    }
	  case shortestPath::SSS :
	    {
	      return std::make_pair(3,0);
	      break;
	    }
	  case shortestPath::SSN :
	    {
	      return std::make_pair(2,1);
	      break;
	    }
	  case shortestPath::SNN :
	    {
	      return std::make_pair(1,2);
	      break;
	    }
	  case shortestPath::NNN :
	    {
	      return std::make_pair(0,3);
	      break;
	    }
	  case shortestPath::NONE :
	    {
	      return std::make_pair(0,0);
	      break;
	    }
	  case shortestPath::AMBIG :
	    {
	      return std::make_pair(SEQMAXUNSIGNED,SEQMAXUNSIGNED);
	      break;
	    }
	  }
	return std::make_pair(SEQMAXUNSIGNED,SEQMAXUNSIGNED);
      }
    catch (SeqException &e)
      {
	throw;
      }
  }

  std::pair<unsigned,shortestPath::pathType> diffType(const std::string &codon1,
						      const std::string &codon2,
						      const Sequence::GeneticCodes & code)
    

  /*!
    \param codon1 a std::string of length 3 representing a codon
    \param codon2 a std::string of length 3 representing a codon
    \param code the genetic code to use in translating the codons
    \return a std::pair<unsigned,shortestPath::pathType>.  The first member of the
    pair takes a value of either 0,1, or 2, depending on the site at which the two codons
    differ (1st, 2nd, or 3rd position, respectively).  If the codons differ 
    at more than 1 site, or contain characters other that
    {A,G,C,T}, the first member will be set to Sequence::SEQMAXUNSIGNED.
    The second member will have the value Sequence::shortestPath::pathType::N
    if the change is nonsynonymous, Sequence::shortestPath::pathType::S if synonymous,
    Sequence::shortestPath::pathType::NONE if the codons don't differ,
    and Sequence::shortestPath::pathType::AMBIG if any of the codons contain characters
    other than {A,G,C,T}.
    \pre (codon1.length()==3 && codon2.length() == 3)
  */
  {
    try
      {
	shortestPath sp(codon1,codon2,code);
	shortestPath::pathType type = sp.type();
	//check preconditions
	if (type == shortestPath::AMBIG)
	  {
	    return std::make_pair(SEQMAXUNSIGNED,type);
	  }
	unsigned site=0,ndiffs=0;
	for(unsigned i = 0 ; i < 3 ; ++i)
	  {
	    if (codon1[i] != codon2[i])
	      {
		++ndiffs;
		site = i;
	      }
	  }
	if( ndiffs != 1 )
	  {
	    return std::make_pair(SEQMAXUNSIGNED,type);
	  }
	return std::make_pair(site,type);
      }
    catch (Sequence::SeqException &e)
      {
	throw;
      }
  }

  boost::tuple<shortestPath::pathType,shortestPath::pathType,shortestPath::pathType>
  diffTypeMulti(const std::string &codon1,
		const std::string &codon2,
		const Sequence::GeneticCodes &code)
    

  /*!
    \return a tuple representing the type of single position changes between codon1
    and codon2.  There is one value in the tuple for each codon position.
    \note The values are assigned as follows:  For each position in codon1, and 2,
    swap the i-th state between the two codons.  If this results in a replacement
    change in both cases, record shortestPath::N.  If it's synonymous in both cases, record
    shortestPath::S.  If the swap results in no change at all (i.e. the two bases are identical),
    record shortestPath::NONE. For all other cases, record shortestPath::AMBIG.  This function
    is most useful at identifying mutations that can be unambiguously classifies as silent
    or replacement.  Note that, if one considers the pathways possible between codons,
    all sites can be assigned as N or S.  For such applications, use Sequence::shortestPath.
   */
  {
    shortestPath::pathType p[3]; //one for each codon position
    std::string t1(codon1),t2(codon2);
    try
      {
	for(unsigned i = 0; i < 3 ; ++i)
	  {
	    //make codons that only differ at one position each
	    std::swap(t1[i],t2[i]);
	    std::pair<unsigned,shortestPath::pathType> dt1 = diffType(t1,codon1,code);
	    std::pair<unsigned,shortestPath::pathType> dt2 = diffType(t2,codon2,code);
	    if (dt1 == dt2)
	      p[i] = dt1.second;
	    else
	      p[i] = shortestPath::AMBIG;
	    //swap back...
	    std::swap(t1[i],t2[i]);
	  }
	return boost::make_tuple(p[0],p[1],p[2]);
      }
    catch(Sequence::SeqException &e)
      {
	throw;
      }
  }
}
