/*!   \include codons.cc
*/

//C++ standard library headers
#include<string>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<functional>
#include<utility>
#include<vector>

//libsequence headers
#include <Sequence/Fasta.hpp>
#include <Sequence/CodonTable.hpp>
#include <Sequence/CountingOperators.hpp>
using namespace std;
using namespace Sequence;

typedef string::size_type sst;

//predicate function to count GC content
struct GorC : public unary_function< char, bool >
{
  inline bool operator()(const char &ch) const
  {
    return ( ch == 'G' || ch == 'C' );
  }
};

int main(int argc, char *argv[])
{
  Fasta sequence;
  CodonUsageTable totalCodonUsage;
  while (! cin.fail())
    {
      cin >> sequence;
      unsigned length = sequence.length();

      //Sequence::Seq has standard iterators pointing to the beginning
      //and end of a sequence.  The functions are Seq::begin() and Seq::end().
      //The iterators are of type std::string::iterator 
      //(and std::string::const_iterator is provided for const Sequence::Seq objects)

      //count the # of gap characters in the sequence
      unsigned numgaps = count(sequence.begin(),sequence.end(),'-');
      //count the # of missing characters in the sequence
      unsigned num_missing = count(sequence.begin(),sequence.end(),'N');

      unsigned datalen = length-numgaps-num_missing;

      //calculate GC content for the whole gene
      unsigned G_or_C = 0;
      G_or_C = count_if(sequence.begin(),sequence.end(),GorC());
      double GC = double(G_or_C)/double(datalen);
	
      //count codon usage
      CodonUsageTable UsageTable= makeCodonUsageTable(&sequence);

      //Use a template operator += from Sequence/CountingOperators.hpp
      //to keep a running total of codon usage
      totalCodonUsage += UsageTable;
      //The class Sequence::Seq has a member function GetSeq(void)
      //which returns the sequence as a std::string.  This allows
      //one to use the rich set of C++ functions for finding 
      //substrings (string::find(), string::rfind(), etc)
      //without having to write wrappers for these functions in the 
      //Seq class
      string seq = sequence.GetSeq();

      //use std::string functions to count GC content at third positions (GC3)
      sst pos = 0;
      unsigned G3=0,C3=0;
      while ( (pos = seq.find("G",pos)) != string::npos)
	{
	  //need to check against pos + 1 because pos is an index,
	  //so adding 1 has the effect of counting positions starting
	  //from 1, so that 3rd positions are 3,6,9,... instead
	  //of 2,5,8,...
	  if (unsigned(pos+1) % 3 == 0.)
	    {
	      ++G3;
	    }
	  ++pos;
	}

      pos = 0;
      while ( (pos = seq.find("C",pos)) != string::npos )
	{
	  if (unsigned(pos+1) % 3 == 0.)
	    {
	      ++C3;
	    }
	  ++pos;
	}

      //divide sum by # of third positions
      double GC3 = double(G3+C3)/( double(datalen) / 3. );

      //output
      cout << "Codon Usage Table for: "<< sequence.GetName() << endl;
      cout << "Sequence length is: "<< sequence.length()<<endl;
      cout << "Length exluding gap characters and missing data is: "<<datalen<<endl;
      cout << "GC content of coding sequence is: "<<GC<<endl;
      cout << "GC3 content of coding sequence is: "<<GC3<<endl;
      unsigned pretty = 1;
      for(unsigned i = 0 ; i < UsageTable.size() ; ++i)
	{
	  if (pretty < 8)
	    {
	      cout << UsageTable[i].first << '\t' << UsageTable[i].second << '\t';
	    }
	  else
	    {
	      cout << UsageTable[i].first << '\t' << UsageTable[i].second << endl;
	      pretty = 0;
	    }
	  ++pretty;
	}
      cout << "//"<<endl;
    }
  cout << "Codon Usage Table for all regions\n";
  unsigned pretty = 1;
  for(unsigned i = 0 ; i < totalCodonUsage.size() ; ++i)
    {
      if (pretty < 8)
	{
	  cout << totalCodonUsage[i].first << '\t' << totalCodonUsage[i].second << '\t';
	}
      else
	{
	  cout << totalCodonUsage[i].first << '\t' << totalCodonUsage[i].second << endl;
	  pretty = 0;
	}
      ++pretty;
    }
  cout << "//\n";
}
