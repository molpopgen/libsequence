#include <Sequence/Fasta.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/SeqUtilities.hpp>
#include <algorithm>
#include <iostream>
#include <functional>

using namespace std;
using namespace Sequence;

/*! \include baseComp.cc */

//an example of using iterators for objects of type Sequence::Seq in conjunction with
//STL algorithms

//usage: cat file | baseComp

struct outputFreq : public std::binary_function<std::pair<char,unsigned>,unsigned,void>
{
  inline void operator ()(const std::pair<char,unsigned> &p, const unsigned &u) const
  {
    std::cout << p.first << " (" << double(p.second)/double(u) << ") ";
  }
};

int main ()
{
  Fasta x;
  while(!cin.fail())
    {
      cin >> x;

      //method 1--use Sequence::stateCounter
      //count base composition for the sequence using the STL algorithm for_each
      //note the use of the iterators belonging to class Sequence::Seq
      //class Sequence::stateCounter is used to keep track of the counts of each base
      stateCounter count = for_each(x.begin(),x.end(),stateCounter('-'));
      cout << "Base composition for sequence "<<x.GetName()<<"\n";

      if (count.ndna == false)
	{
	  unsigned len = x.length() - count.gap - count.n;//ungapped length,exclude missing data
	  
	  //turn the counts into percentages
	  double percA = (count.a>0) ? double(count.a)/double(len) : 0.;
	  double percG = (count.g>0) ? double(count.g)/double(len) : 0.;
	  double percC = (count.c>0) ? double(count.c)/double(len) : 0.;
	  double percT = (count.t>0) ? double(count.t)/double(len) : 0.;
	  
	  //output
	  std::cout << "using Sequence::stateCounter: "
		    << "A (" << percA << ") "
		    << "G (" << percG << ") "
		    << "C (" << percC << ") "
		    << "T (" << percT << ")\n";

	}
      else
	{
	  cout << "non-DNA character encountered.  Skipping...\n";
	}
      //method 2--use Sequence::makeCountList (doesn't check for non-standard DNA characters)
      std::map<char,unsigned> m = Sequence::makeCountList(x.begin(),x.end());
      std::cout << "using Sequence::makeCountList(): ";
      std::for_each(m.begin(),m.end(),std::bind2nd(outputFreq(),x.length()));
      std::cout << std::endl;
    }
}
