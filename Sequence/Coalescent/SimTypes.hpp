#ifndef __SEQUENCE_COALESCENT_SIM_TYPES_HPP__
#define __SEQUENCE_COALESCENT_SIM_TYPES_HPP__

#include <vector>
#include <list>
#include <iosfwd>
#include <cassert>
//#include <memory>
/*! \file SimTypes.hpp
  @brief declaration of types for coalescent simulation
 */

/*! \struct Sequence::segment Sequence/Coalescent/SimTypes.hpp
  @brief A portion of a recombining chromosome

  A segment is a portion of a recombining chromosome. A segment
  is defined in terms of its mutational length.  Therefore,
  a segment consisting of k sites has k-1 positions at
  which recombination events can occur.
  
  This struct contains 3 public members: 
  - beg:  the first site in the segment
  - end:  the last site in the segment
  - desc: the individual in the sample to which the segment leads
  
  \ingroup coalescent
*/

/*! \struct Sequence::chromosome Sequence/Coalescent/SimTypes.hpp
  @brief A chromosome is a container of segments.
  \note this is a malloc-based container
  \ingroup coalescent
*/

/*! \struct Sequence::node Sequence/Coalescent/SimTypes.hpp
  @brief A point on a marginal tree at which a coalescent event occurs.

  A node is a branch-point on a coalescent tree.  This struct
  contains two public members.
  \ingroup coalescent
*/

/*! \struct Sequence::marginal Sequence/Coalescent/SimTypes.hpp
  @brief The genealogy of a portion of a chromosome on which no
  recombination has occurred.

  A marginal history is a coalscent tree for a region
  in which no recombination has occured in the history of
  a sample.  A sorted list of marginal trees is the 
  ancstral recombination graph for a recombining region,
  and the typedef std::list<marginal> arg is used in this
  library
  \ingroup coalescent
*/

/*! \class Sequence::newick_stream_marginal_tree Sequence/Coalescent/SimTypes.hpp
  @brief Class that provides a typecast-on-output of a marginal tree to a newick tree
  Example use:
  \code
  //assume we've done something useful to get
  //data into an arg called sample_history
  //tlength is the total length of the region being simulated
  int seg = 0,nsegs=sample_history.size();
  arg::iterator i = sample_history.begin(),j=i;
  ++j;
  for( seg=0 ; seg < nsegs ; ++seg,++i,++j)
  {
  int length = ( (seg<nsegs-1) ? (j->beg-i->beg) : (tlength-i->beg) );
  cout << '[' 
  << length
  << ']'
  << newick_stream_marginal_tree(i)
  << '\n';
  }
  \endcode
  \ingroup coalescent
*/

namespace Sequence
{
  struct segment
  {
    int beg,end,desc;
    segment();
    explicit segment(const int & b,
		     const int & e,
		     const int & d);
  };

  struct chromosome
  {
    /*!
      The list of segments making up the ancestral material
      of a chromosome at the current point in the simulation.
      \note Please be aware that allocation of this array
      is done via malloc, rather than operator new.  As coalescent 
      simulations are all about copying segments, the use of malloc is
      more efficient for simulation purposes, as it results
      in many fewer calls to the default constructor.  But
      this means you cannot assign to this pointer something allocated 
      with operator new, else bad things will happen.
    */
    segment * segs;
    typedef segment * iterator;
    typedef const segment * const_iterator;
    /*!
      the population in which the chromosome is currently found
    */
    int pop;
    /*!
      The number of segments contained in the pointer segs
    */
    unsigned nsegs;
    chromosome();
    chromosome(const chromosome & ch);
    chromosome( const std::vector<segment> & initial_segs,
		const int & population = 0 );
    ~chromosome();
    chromosome & operator=(const chromosome & ch);
    void swap_with( chromosome & ch );
    void assign_allocated_segs( segment * newsegs,
				const int & new_nsegs );
    int first() const
    /*!
      \return the first position in the chromosome
    */
    {
      assert(nsegs>0);
      return segs->beg;
    }

    int last() const
    /*!
      \return the last position in the chromosome
    */
    {
      assert(nsegs>0);
      return (segs+nsegs-1)->end;
    }
    int links() const;
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
  };

  std::ostream & operator<<(std::ostream & s,const  chromosome & c);

  struct node
  {
    /*!
      The (coalescent-scaled) time at which the node was generated
    */ 
    double time;
    /*!
      The index in the marginal tree that is the ancestor 
      of the current node
    */
    int abv;
    node(const double & t = 0.,
	 const int & a = -1 );
  };

  struct marginal
  {
    /*! The (mutational) site at which the current
      marginal tree begins
    */
    int beg;
    /*!
      The sample size being simulated.  The 2*nsam-1 nodes
      in the tree are therefore indexed from 0 to 2*nsam-2
    */
    mutable int nsam;
    /*!
      The current number of nodes in the tree.  By current, it 
      is meant the current time in the simulation.  At the start
      of a simulation of a sample size of k chromosomes, this should
      be initialized to k.  This value is manipulated as nodes
      are added into the marginal history by the coalesce function.
    */
    int nnodes;
    /*!
      tree is the coalescent history of this marginal tree
    */
    std::vector<node> tree;
    typedef std::vector<node>::iterator iterator;
    typedef std::vector<node>::const_iterator const_iterator;
    typedef std::vector<node>::size_type size_type;
    typedef std::vector<node>::reference reference;
    typedef std::vector<node>::const_reference const_reference;
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    reference operator[](const std::vector<node>::size_type &i);
    const_reference operator[](const std::vector<node>::size_type &i) const;
    marginal(const int & b, const int & ns,const int & nn, 
	     const std::vector<node> & tree);
    bool operator<(const marginal & rhs) const;
  };

  /*!
    @brief Ancestral Recombination Graph

    An arg is an "ancestral recombination graph",
    which is a linked list of marginal histories.

    \note The implementation of the crossover function
    ensures that the marginal trees are sorted in ascending
    order determined by marginal::beg

    \ingroup coalescent
  */
  typedef std::list<marginal> arg;
  std::ostream & operator<<(std::ostream & s, const marginal & m);

  class newick_stream_marginal_tree_impl;
  class newick_stream_marginal_tree
  {
  private:
    newick_stream_marginal_tree_impl * impl;
    /*
    marginal::const_iterator mi;
    const int nsam;
    std::vector<int> left,right;
    std::vector<node> tree;
    void init();
    std::ostream & parens( const int & noden,
			   std::ostream & o) const;
    */
  public:
    newick_stream_marginal_tree( const marginal & m );
    newick_stream_marginal_tree( const marginal * m );
    newick_stream_marginal_tree( arg::const_iterator m );
    newick_stream_marginal_tree( arg::iterator m );
    ~newick_stream_marginal_tree( );
    std::vector<node> get_tree() const;
    std::ostream & print( std::ostream & o ) const;
    std::istream & read( std::istream & i );
  };
  std::ostream & operator<<(std::ostream & o, const newick_stream_marginal_tree & n);
  std::istream & operator>>(std::istream & i, newick_stream_marginal_tree & n);
}
#endif
