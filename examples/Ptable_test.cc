/*! \include Ptable_test.cc
 */

#include <Sequence/Ptable.hpp>
#include <Sequence/PolySites.hpp>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Sequence;

int main( int argc, char ** argv )
{
  //Can use C++11 initializer lists
  //polymorphicSite is a typedef for pair<double,string>
  Ptable x{ polymorphicSite(1.,"AAAT"),polymorphicSite(2.,"GAGG") };

  //Print the data
  cout << "The data:\n";
  for( auto _x : x )
    {
      cout << _x.first << ' ' << _x.second << '\n';
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
  
  //Delete sites (this means that you can use algorithms + function objects on all polymorphic sites!)
  x.erase(x.begin());

  cout << "The data after erasing:\n";
  for( auto _x : x )
    {
      cout << _x.first << ' ' << _x.second << '\n';
    }

  //Convert to a PolySites
  PolySites ps(x.begin(),x.end());

  cout << "Print out the classic object:\n"
       << ps << '\n';
}
