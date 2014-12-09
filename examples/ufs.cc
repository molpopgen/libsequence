/*
  ufs.cc
  
  Calculate the polarized (unfolded) site frequency spectrum from a file containing
  aligned sequences in Fasta format.

  Author: Kevin Thornton
*/

#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SeqUtilities.hpp>
#include <Sequence/CountingOperators.hpp>
#include <algorithm>
#include <vector>
#include <iostream>
#include <functional>
#include <cctype>

//so that you know where things come from
using std::map;
using std::vector;
using std::cout;
using std::cerr;
using Sequence::Fasta;
using Sequence::PolySites;
using Sequence::makeCountList; // <Sequence/SeqUtilities.hpp>
using Sequence::operator+=;    // <Sequence/CountingOperators.hpp> defines an operator+= for maps
using Sequence::Alignment::GetData;
using Sequence::Alignment::IsAlignment;

bool validStates(const std::map<char,unsigned> & counts)
/*
  counts is assumed to be a map of nucleotides and their number of occurences.
  The function returns false if any of the nucleotides are not in the set {A,G,C,T,N},
  and toupper is used so that matching is case-insensitive.
 */
{
  for( map<char,unsigned>::const_iterator i = counts.begin() ;
       i != counts.end() ;
       ++i )
    {
      char ch = toupper(i->first);
      if ( ch != 'A' && ch != 'G' && ch != 'C' && ch != 'T' && ch != 'N' )
	return false;
    }
  return true;
}

int main(int argc, char **argv)
{
  if(argc!=3)
    {
      cerr << "usage: ufs file outgroup_index\n";
      exit(1);
    }

  const char * fastafile = argv[1];
  const unsigned outgroup = atoi(argv[2]);

  //read in alignment and make sure that all elements are the same length, otherwise exit
  vector<Fasta> alignment;
  GetData(alignment,fastafile);
  if( ! IsAlignment(alignment) )
    {
      cerr << fastafile << " not aligned\n";
      exit(1);
    }

  //make a table of variable sites, excluding gaps and sites with > 2 states in the alignment
  PolySites SNPtable(alignment,true); //the "true" eliminates all sites in alignment with > 2 states

  vector< unsigned > ufs(alignment.size()-1,0); //store the frequency spectrum

  //here, we use iterators to polymorphic sites to access columns in the SNP table
  //these iterators are typedefs for std::pair<double,std::string>
  for( PolySites::const_site_iterator i = SNPtable.sbegin() ;
       i != SNPtable.send() ; ++i )
    {
      const char ancstate = i->second[outgroup];

      //now, we want to count up all the nucelotide states in the ingroup.  we do this in 2 stages,
      //so that we guarantee that we skip the outgroup.  as the second element in this iterator
      //is a std::string, normal STL range operations apply.

      //the function Sequence::makeCountList takes two iterators, of type ITR, (beg,end) as arguments,
      //and counts all the elements in the range beg,end-1, returning a 
      //std::map< std::iterator_traits<ITR>::value_type, unsigned >.  In this case, the 
      //value_type is a char, since the iterators are those of std::string
      std::map<char,unsigned> counts = makeCountList(i->second.begin(),i->second.begin()+outgroup);

      //here is a trick.  std::map does not have an operator+=, but <Sequence/CountingOperators.hpp>
      //provides template to define these operators, which is useful for counting operations like this
      counts += makeCountList(i->second.begin()+outgroup+1,i->second.end());

      if( ! validStates(counts) )
	{
	  cerr << "site " << i->first << " contains characters other than A,G,C,T,N\n";
	}
      else
	{
	  for( map<char,unsigned>::const_iterator i = counts.begin() ;
	       i != counts.end() ;
	       ++i )
	    {
	      if ( toupper(i->first) != ancstate )
		{
		  ufs[i->second-1]++;
		}
	    }
	}
    }

  //output the unfolded site frequency spectrum.
  for(unsigned i=0;i<ufs.size();++i)
    {
      cout << i+1 << '\t' << ufs[i] << '\n';
    }
  return 0;
}
