/*! \file TreeOperations.hpp
  Things you may want to do with marginal trees in a coalescent simulations
 */
#ifndef __SEQUENCE_COALESCENT_TREE_OPERATIONS_HPP__
#define __SEQUENCE_COALESCENT_TREE_OPERATIONS_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <vector>
#include <memory>
namespace Sequence
{
  double total_time( const marginal::const_iterator beg,
		     const int & nsam );
  
  int pick_branch( marginal::const_iterator beg,
		   const int & nsam,
		   const double & rtime);
  
  std::vector<int> get_all_descendants (marginal::const_iterator beg,
					const int & nsam,
					const int & branch);

  bool is_descendant( marginal::const_iterator beg,
		      const int & ind,
		      const int & branch );

  double total_time_on_arg( const Sequence::arg & sample_history,
			    const int & total_number_of_sites );

  void minimize_arg( arg * sample_history );

  class sfs_times_impl;
  class sfs_times
  {
  private:
    std::auto_ptr<sfs_times_impl> impl;
  public:
    sfs_times();
    sfs_times(arg::const_iterator sample_history_beg,
	      const arg::size_type & nsegs,
	      const int & total_nsites_simulated,
	      bool folded = false);
    sfs_times(const sfs_times &);
    ~sfs_times();
    
    double operator[]( std::vector<double>::size_type const & ) const;
    sfs_times & operator=(const sfs_times &);
    bool operator==(const sfs_times & rhs) const;
    double ttime() const;
    size_t size() const;
    typedef std::vector<double>::const_iterator const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
  };
}
#endif
