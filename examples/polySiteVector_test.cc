/*! \include polySiteVector_test.cc
  \file polySiteVector_test.cc @brief Examples of using Sequence::polySiteVector
 */

#include <Sequence/polySiteVector.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/stateCounter.hpp>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;
using namespace Sequence;

//uses the "classic" stateCounter to determine MAF
double local_maf( const polymorphicSite & p )
{
  stateCounter sc = for_each(p.second.begin(),
			     p.second.end(),stateCounter());
  unsigned mcount = numeric_limits<unsigned>::max();
  mcount = sc.a ? min(mcount,sc.a) : mcount;
  mcount = sc.g ? min(mcount,sc.g) : mcount;
  mcount = sc.c ? min(mcount,sc.c) : mcount;
  mcount = sc.t ? min(mcount,sc.t) : mcount;
  return double(mcount)/double(p.second.size()-sc.n);
}

std::vector<std::string::size_type> local_has_missing( const polymorphicSite & p )
//case-insensitive!
{
  typedef std::vector<std::string::size_type> rvtype;
  rvtype rv;

  string::size_type f = p.second.find('N');
  while (f != string::npos)
    {
      rv.push_back(f);
      f = p.second.find('N',f+1);
    }
  return rv;
}

int main( int argc, char ** argv )
{
  //Can use C++11 initializer lists
  //polymorphicSite is a typedef for pair<double,string>
  polySiteVector x{ polymorphicSite(1.,"AAAT"),polymorphicSite(2.,"GANG") };

  //Print the data
  cout << "The data:\n";
  for( auto _x : x )
    {
      cout << _x.first << ' ' << _x.second << '\n';
    }

  //Calculate MAFs
  vector<double> mafs;
  for_each( x.begin(),x.end(), [&mafs](const polymorphicSite & __p){ mafs.push_back(local_maf(__p)); } );
  cout << "The minor allele frequencies:\n";
  for( auto __maf : mafs )
    {
      cout << __maf << '\n';
    }

  //Reverse the order of the segregating sites by sort.
  //With "classic" PolyTables, this could not be done, b/c access to polymorphicSites was const-only
  sort(x.begin(),x.end(),[](const polymorphicSite & lhs,
			    const polymorphicSite & rhs) { return lhs.first > rhs.first; });

  cout << "The data after sorting:\n";
  for( auto _x : x )
    {
      cout << _x.first << ' ' << _x.second << '\n';
    }

  //Who has missing data?
  std::vector<std::string::size_type> missing;
  for_each(x.begin(),x.end(),[&missing](const polymorphicSite & __p) { 
      std::vector<std::string::size_type> __temp = local_has_missing(__p);
      copy(__temp.begin(),__temp.end(),back_inserter(missing));
    });

  cout << "Index of individuals with missing data:\n";
  for( auto __m : missing )
    {
      cout << __m << '\n';
    }


  //Delete sites with MAF <= 0.25
  x.erase( remove_if(x.begin(),x.end(),[](const polymorphicSite & __p){ return local_maf(__p) <= 0.25; } ),x.end() );

  cout << "The data after applying MAF filter:\n";
  for( auto _x : x )
    {
      cout << _x.first << ' ' << _x.second << '\n';
    }

  //Convert to a PolySites
  PolySites ps(x.begin(),x.end());

  cout << "Print out the classic object:\n"
       << ps << '\n';
}
