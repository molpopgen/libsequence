//  -*- C++ -*- 
#ifndef __SEQUENCE_COALESCENT_BITS_RECOMBINATION_TCC__
#define __SEQUENCE_COALESCENT_BITS_RECOMBINATION_TCC__
#include <utility>
#include <numeric>

namespace Sequence
{
#ifndef DOXYGEN_SKIP
  template<typename uniform01_generator>
  std::pair<int,int> pick_spot_details(uniform01_generator & uni01,
				       const double & total_reclen,
				       const std::vector<double> & reclens,
				       std::vector<chromosome>::const_iterator sample_begin,
				       const unsigned & current_nsam,
				       const double * rec_map)
  {
    double sum=0.;
    int recombinant=0;
    int pos=0;
    double ran = uni01();
    bool flag = false;
    while( (sample_begin+recombinant) < 
	   (sample_begin+current_nsam) && ! flag)
      {
	sum += reclens[recombinant];
	if ( ran <= sum/total_reclen )
	  {
	    int beg = (sample_begin+recombinant)->segs->beg;
	    int end = ((sample_begin+recombinant)->segs+(sample_begin+recombinant)->nsegs-1)->end;
	    double cum_prob = std::accumulate(rec_map+beg,rec_map+end,0.);
	    double ran2 = uni01(),sum2=0.;
	    //we only want to simulate crossover spots
	    //from the density that covers the interval beg to end-1
	    for(int i=beg ; i < end ; ++i)
	      {
		sum2+= *(rec_map+i);
		if ( ran2 <= sum2/cum_prob )
		  {
		    pos = i;
		    flag = true;
		    break;
		  }
	      }
	  }
	if(flag) break;
	++recombinant;
      }
    return std::make_pair(recombinant,pos);
  }
#endif

  template<typename uniform01_generator>
  std::pair<int,int> pick_spot( uniform01_generator & uni01,
				const double & total_reclen,
				const std::vector<double> & reclens,
				std::vector<chromosome>::const_iterator sample_begin,
				const unsigned & current_nsam,
				const double * rec_map)
  /*!
    Picks a positions  amongst all chromosomes at which a recombination event
    will occur, based on an arbitrary genetic map
    \param uni01 a function/object which takes no arguments and can return a U[0,1)
    \param total_reclen the total recombination length of all chromosomes in the sample
    \param reclens a vector of the proportion of \a total_reclen contributed by each chromosome.
    This needs to be ordered in the same order as \a sample_begin to (\a sample_begin + \a current_nsam - 1)
    \param sample_begin an iterator pointing to the beginning of the sample
    \param current_nsam the current sample size in the simulation
    \param rec_map an array of probabilities describing the recombination map.  The map is completely up to
    the programmer, and it is not checked for sanity at all in this function.  For a region of k sites,
    indexes 0 to k-2 of this array should be filled.  The i-th element should contain the probability
    that a crossover occurs between position i and i+1.  The sum of all elements should be 1, such 
    that the array describes the recombination map in terms of a probability distribution function. An 
    example of how to do this is in the file examples/msbeta.cc that comes with the source for this library.
    \return a pair of integers containing the index of the recombinant chromosome (.first),
    and the position at which the crossover will occur (.second)
    \ingroup coalescent
  */
  {
    return pick_spot_details(uni01,total_reclen,reclens,sample_begin,current_nsam,rec_map);
  }

  template<typename uniform01_generator>
  std::pair<int,int> pick_spot( const uniform01_generator & uni01,
				const double & total_reclen,
				const std::vector<double> & reclens,
				std::vector<chromosome>::const_iterator sample_begin,
				const unsigned & current_nsam,
				const double * rec_map)
  /*!
    Picks a positions  amongst all chromosomes at which a recombination event
    will occur, based on an arbitrary genetic map
    \param uni01 a function/object which takes no arguments and can return a U[0,1)
    \param total_reclen the total recombination length of all chromosomes in the sample
    \param reclens a vector of the proportion of \a total_reclen contributed by each chromosome.
    This needs to be ordered in the same order as \a sample_begin to (\a sample_begin + \a current_nsam - 1)
    \param sample_begin an iterator pointing to the beginning of the sample
    \param current_nsam the current sample size in the simulation
    \param rec_map an array of probabilities describing the recombination map.  The map is completely up to
    the programmer, and it is not checked for sanity at all in this function.  For a region of k sites,
    indexes 0 to k-2 of this array should be filled.  The i-th element should contain the probability
    that a crossover occurs between position i and i+1.  The sum of all elements should be 1, such 
    that the array describes the recombination map in terms of a probability distribution function. An 
    example of how to do this is in the file examples/msbeta.cc that comes with the source for this library.
    \return a pair of integers containing the index of the recombinant chromosome (.first),
    and the position at which the crossover will occur (.second)
    \ingroup coalescent
  */
  {
    return pick_spot_details(uni01,total_reclen,reclens,sample_begin,current_nsam,rec_map);
  }
}//namespace Sequence
#endif
