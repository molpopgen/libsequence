#ifndef __SEQUENCE_COALESCENT_COALESCE_HPP__
#define __SEQUENCE_COALESCENT_COALESCE_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <utility>
#include <boost/ref.hpp>
namespace Sequence
{
  template<typename uniform_generator>
  std::pair<int,int> pick2_in_deme( uniform_generator & uni, 
				    const std::vector<Sequence::chromosome> & sample,
				    const int & ttl_nsam,
				    const int & deme_nsam,
				    const int & deme );

  template<typename uniform_generator>
  std::pair<int,int> pick2_in_deme( const uniform_generator & uni, 
				    const std::vector<Sequence::chromosome> & sample,
				    const int & ttl_nsam,
				    const int & deme_nsam,
				    const int & deme );

  template<typename uniform_generator>
  std::pair<int,int> pick2( uniform_generator & uni, const int & nsam);

  template<typename uniform_generator>
  std::pair<int,int> pick2( const uniform_generator & uni, const int & nsam);

  bool isseg( chromosome::const_iterator seg, const int nsegs,
	      const int pos, int * offset );

  int coalesce(const double & time,
	       const int & ttl_nsam,
	       const int & current_nsam,
	       const int & c1,
	       const int & c2,
	       const int & nsites,
	       int * nlinks,
	       std::vector<chromosome> * sample,
	       arg * sample_history);
}
#endif
#include <Sequence/Coalescent/bits/Coalesce.tcc>
