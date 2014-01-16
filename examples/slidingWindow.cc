#include <Sequence/PolySites.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <vector>
#include <iostream>

/*! \include slidingWindow.cc */

//Read in a data set of aligned sequence in Fasta
//format.  Create a polymorphism table.  Calculate
//Tajima's D for the whole table.  Then, run a sliding
//window of 1 segregating site (with a jump size of 1)
//along the SNP table, and use that to calculate Tajima's
//D for each site.

//This is a somewhat contrived example, but it illustrates
//the sliding window code.

int main(int argc, char **argv)
{
  const char * infilename = argv[1];
  std::vector<Sequence::Fasta> data;
  Sequence::Alignment::GetData(data,infilename);

  if ( Sequence::Alignment::IsAlignment(data) && 
       Sequence::Alignment::validForPolyAnalysis(data.begin(),data.end()) )
    {
      Sequence::PolySites SNPtable(data);
      if (! SNPtable.empty())
	{
	  Sequence::PolySNP analyzeRegion(&SNPtable);
	  std::cout << "Tajima's D for the region is: "<< analyzeRegion.TajimasD() << std::endl;
	  
	  Sequence::PolyTableSlice<Sequence::PolySites> windows(SNPtable.sbegin(),
								SNPtable.send(),1,1);
	  Sequence::PolyTableSlice<Sequence::PolySites>::const_iterator itr = windows.begin();
	  while(itr < windows.end())
	    {
	      Sequence::PolySites window = windows.get_slice(itr);
	      Sequence::PolySNP analyzeWindow(&window);
	      std::cout << "D for window " 
			<< itr-windows.begin()
			<< " is: "
			<< analyzeWindow.TajimasD()
			<< std::endl;
	      ++itr;
	    }
	}
    }
}
