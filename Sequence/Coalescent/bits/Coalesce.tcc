//  -*- C++ -*- 
#ifndef __SEQUENCE_COALESCENT_BITS_COALESCE_TCC__
#define __SEQUENCE_COALESCENT_BITS_COALESCE_TCC__

namespace Sequence
{
#ifndef DOXYGEN_SKIP
  template<typename uniform_generator>
  std::pair<int,int> pick2_in_deme_details( uniform_generator &  uni, 
					    const std::vector<Sequence::chromosome> & sample,
					    const int & current_nsam,
					    const int & deme_nsam,
					    const int & deme )
  {
    assert( deme_nsam > 0 );
    assert( deme >= 0 );
    assert( deme_nsam <= current_nsam );
    std::pair<int,int> two = pick2(uni,deme_nsam);
    assert(two.first != two.second);
    int c1=0,c2=0;
    for(int i=0,j=0;i<current_nsam;++i)
      {
	if( (sample.begin()+i)->pop == deme && j == two.first )
	  c1 = i;
	else if ( (sample.begin()+i)->pop == deme && j == two.second )
	  c2 = i;
	if( (sample.begin()+i)->pop == deme ) ++j;
      }
    assert(c1!=c2);
    return std::make_pair(std::min(c1,c2),std::max(c1,c2));
  }

  template<typename uniform_generator>
  std::pair<int,int> pick2_details(uniform_generator & uni, const int & nsam)
  {
    int i = int(uni(boost::cref(0.),boost::cref(nsam)));
    int j;
    while( (j = int(uni(boost::cref(0.),boost::cref(nsam)))) == i )
      {
      }
    if(i>j) std::swap(i,j);
    return std::make_pair(i,j);
  }

  template<typename uniform_generator>
  std::pair<int,int> pick2_in_deme( uniform_generator & uni, 
				    const std::vector<Sequence::chromosome> & sample,
				    const int & current_nsam,
				    const int & deme_nsam,
				    const int & deme )
  /*!
    \brief Choose two random chromosomes from the same deme
    \param uni A random number generator taking two arguments, a and b, and returning 
    a random variable distributed uniformly over [a,b)
    \param sample The current state of the simulated sample
    \param current_snam The total sample size being simuled (the sum of sample sizes over all demes)
    \param deme_nsam The sample size of the deme from which you wish to sample
    \param deme The index ( 0 <= deme < # populations ) of the deme from which you wish to sample 
    \return A pair of integers which contains the indexes of two chromosomes in @a sample
    \ingroup coalescent
  */
  {
    return pick2_in_deme_details(uni,sample,current_nsam,deme_nsam,deme);
  }
#endif

  template<typename uniform_generator>
  std::pair<int,int> pick2_in_deme( const uniform_generator & uni, 
				    const std::vector<Sequence::chromosome> & sample,
				    const int & current_nsam,
				    const int & deme_nsam,
				    const int & deme )
  /*!
    \brief Choose two random chromosomes from the same deme
    \param uni A random number generator taking two arguments, a and b, and returning 
    a random variable distributed uniformly over [a,b)
    \param sample The current state of the simulated sample
    \param current_nsam The total sample size being simuled (the sum of sample sizes over all demes)
    \param deme_nsam The sample size of the deme from which you wish to sample
    \param deme The index ( 0 <= deme < # populations ) of the deme from which you wish to sample 
    \return A pair of integers which contains the indexes of two chromosomes in @a sample
    \ingroup coalescent
  */
  {
    return pick2_in_deme_details(uni,sample,current_nsam,deme_nsam,deme);
  }

  template<typename uniform_generator>
  std::pair<int,int> pick2( uniform_generator & uni, const int & nsam)
  /*!
    \param uni a random number function/object capable of returning
    a double-precision random number between 0 and \a nsam-1
    \param nsam the current sample size in the simulation
    \return A pair of integers which contains the indexes of two chromosomes in @a sample
    \ingroup coalescent
  */
  {
    return pick2_details(uni,nsam);
  }

  template<typename uniform_generator>
  std::pair<int,int> pick2( const uniform_generator & uni, const int & nsam)
  /*!
    \param uni a random number function/object capable of returning
    a double-precision random number between 0 and \a nsam-1
    \param nsam the current sample size in the simulation
    \return A pair of integers which contains the indexes of two chromosomes in @a sample
    \ingroup coalescent
  */
  {
    return pick2_details(uni,nsam);
  }
}

#endif
