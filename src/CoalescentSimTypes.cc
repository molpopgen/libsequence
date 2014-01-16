/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <Sequence/Coalescent/SimTypes.hpp>
#include <cassert>
#include <iostream>
#include <cstdlib>

namespace Sequence
{
  segment::segment() : beg(0),end(0),desc(0)
		  /*!
		    @brief default constructor
		    initializes all members to 0
		  */
  {
  }

  segment::segment(const int & b,
		   const int & e,
		   const int & d = 0) : beg(b),end(e),desc(d)
		  /*!
		    @brief constructor
		    \param b beg
		    \param e end
		    \param d desc
		  */
  {
    assert ( end >= beg );
  }

  chromosome::chromosome() : segs(NULL),
			     pop(0),
			     nsegs(0)
			/*! 
			  @brief constructor
			  sets segs to NULL, pop to 0, and nsegs to 0
			*/
  {
  }

  chromosome::chromosome( const std::vector<segment> & initial_segs,
			  const int & population ) : 
    segs((segment *)malloc(initial_segs.size()*
			   sizeof(segment))), 
    pop(population),
    nsegs(int(initial_segs.size()))
			/*!
			  @brief constructor
			  \param initial_segs a vector of segments
			  \param population used to set pop
			  \note nsegs is set to initial_segs.size()
			*/
  {
    std::copy(initial_segs.begin(),initial_segs.end(),segs);
  }

  chromosome::~chromosome()
  /*!
    frees pointer to segments
  */
  {
    if(nsegs>0) free(segs);
  }

  chromosome::chromosome( const chromosome & ch ) :
    segs((segment *)malloc(ch.nsegs*sizeof(segment))),
    pop(ch.pop),
    nsegs(ch.nsegs)
			/*!
			  @brief copy constructor
			*/
  {
    std::copy(ch.segs,ch.segs+ch.nsegs,segs);
  }

  chromosome & chromosome::operator=(const chromosome & ch)
  /*!
    @brief assignment operator
  */
  {
    if(this == &ch) return *this;
    if(ch.nsegs>this->nsegs)
      {
	this->segs = (segment *)realloc(this->segs,ch.nsegs*sizeof(segment));
      }
    std::copy(ch.segs,ch.segs+ch.nsegs,segs);
    this->nsegs = ch.nsegs;
    this->pop = ch.pop;
    return *this;
  }


  void chromosome::swap_with( chromosome & ch )
  /*!
    Swaps the data members of the current chromosome with chromosome ch.
    Called by the coalesce routine, and is necessary to prevent nastiness
    such as multiple calls to free when vectors of chromosomes go out of scope.
    Implemented as:
    std::swap(this->segs,ch.segs);
    std::swap(this->nsegs,ch.nsegs);
    std::swap(this->pop,ch.pop);
  */
  {
    std::swap(this->segs,ch.segs);
    std::swap(this->nsegs,ch.nsegs);
    std::swap(this->pop,ch.pop);
  }

  void chromosome::assign_allocated_segs( segment * newsegs,
					  const int & new_nsegs )
  /*!
    Replaces the current segs with those pointed to by newsegs
    \param newsegs an array of segments allocated with malloc
    \param new_nsegs the number of segs stored in \a newsegs
  */
  {
    free(segs);
    segs=newsegs;
    nsegs=new_nsegs;
  }

  chromosome::iterator chromosome::begin()
  /*!
    \return segs
  */
  {
    return segs;
  }

  chromosome::iterator chromosome::end()
  /*!
    \return segs+nsegs
  */
  {
    return segs+nsegs;
  }

  chromosome::const_iterator chromosome::begin() const
  /*!
    \return segs
  */
  {
    return segs;
  }

  chromosome::const_iterator chromosome::end() const
  /*!
    \return segs+nsegs
  */
  {
    return segs+nsegs;
  }

  int chromosome::links() const
  /*!
    Computes and returns the number of positions
    at which recombination can occur in the chromosome
    \return  (segs+nsegs-1)->end - segs->beg
  */
  {
    return( (segs+nsegs-1)->end - segs->beg );
  }

  node::node(const double & t,
	     const int & a) : time(t),abv(a)
	    /*!
	      \param t time
	      \param a abv
	    */
  {
  }

  marginal::marginal(const int & b, const int & ns,
		     const int & nn,
		     const std::vector<node> & t) :
    beg(b),nsam(ns),nnodes(nn),tree(t)
		    /*!
		      \param b beg
		      \param ns nsam
		      \param nn nnodes
		      \param t tree
		    */
  {
  }

  marginal::iterator marginal::begin()
  /*!
    \return this->tree.begin();
  */
  {
    return tree.begin();
  }

  marginal::iterator marginal::end()
  /*!
    \return this->tree.end();
  */
  {
    return tree.end();
  }

  marginal::const_iterator marginal::begin() const
  /*!
    \return this->tree.begin();
  */
  {
    return tree.begin();
  }

  marginal::const_iterator marginal::end() const
  /*!
    \return this->tree.end();
  */
  {
    return tree.end();
  }

  marginal::reference marginal::operator[](const std::vector<node>::size_type &i)
  /*!
    \return this->tree[i]
    \warning Range checking by assert only!
  */
  {
    assert( int(i) <= 2*nsam-2 );
    return tree[i];
  }

  marginal::const_reference marginal::operator[](const std::vector<node>::size_type &i) const
  /*!
    \return this->tree[i]
    \warning Range checking by assert only!
  */
  {
    assert( int(i) <= 2*nsam-2 );
    return tree[i];
  }

  bool marginal::operator<(const marginal & rhs) const
  /*!
    return this->beg < rhs.beg;
  */
  {
    return this->beg < rhs.beg;
  }

  std::ostream & operator<<(std::ostream & s,const  chromosome & c)
  /*!
    @brief output operator for chromosome types in coalescent simulation
    Outputs the segments contained by the chromosome
    \ingroup operators
  */
  {
    s << '(';
    chromosome::const_iterator seg = c.begin();
    for( ; seg < c.end() ; 
	 ++seg)
      {
	s << seg->beg << ':' << seg->end << ',';
      }
    s << ' ' << c.nsegs << ' ' << c.links() << ')';
    return s;
  }

  std::ostream & operator<<(std::ostream & s, const marginal & m)
  /*!
    @brief Write a marginal tree to an ostream
    \ingroup operators
  */
  {
    s << "marginal tree begins at: " << m.beg << '\n';
    int seg=0;
    for( ;seg < m.nnodes ; ++seg)
      {
	s << seg << ' ' 
	  << (m.begin()+seg)->time << ' '
	  << (m.begin()+seg)->abv << '\n';
      }
    s << seg << ' '
      << (m.begin()+seg)->time << ' '
      << (m.begin()+seg)->abv;
    return s;
  }

  struct newick_stream_marginal_tree_impl
  /*!
    Implementation details of newick_stream_marginal_tree
  */
  {
    marginal::const_iterator mi;
    const int nsam;
    std::vector<int> left,right;
    std::vector<node> tree;
    void init();
    std::ostream & parens( const int & noden,
			   std::ostream & o) const;
    newick_stream_marginal_tree_impl(const marginal & m);
    newick_stream_marginal_tree_impl(const marginal * m);
    newick_stream_marginal_tree_impl(arg::const_iterator m);
    newick_stream_marginal_tree_impl(arg::iterator m);
  };

  newick_stream_marginal_tree_impl::newick_stream_marginal_tree_impl(const marginal & m) :
    mi(m.begin()),nsam(m.nsam),
    left(2*nsam-1,-1),right(left),
    tree(std::vector<node>())
  {
    init();
  }

  newick_stream_marginal_tree_impl::newick_stream_marginal_tree_impl(const marginal * m) :
    mi(m->begin()),nsam(m->nsam),
    left(2*nsam-1,-1),right(left),
    tree(std::vector<node>())
  {
    init();
  }

  newick_stream_marginal_tree_impl::newick_stream_marginal_tree_impl(arg::const_iterator m) :
    mi(m->begin()),nsam(m->nsam),
    left(2*nsam-1,-1),right(left),
    tree(std::vector<node>())
  {
    init();
  }

  newick_stream_marginal_tree_impl::newick_stream_marginal_tree_impl(arg::iterator m) :
    mi(m->begin()),nsam(m->nsam),
    left(2*nsam-1,-1),right(left),
    tree(std::vector<node>())
  {
    init();
  }

  std::ostream &
  newick_stream_marginal_tree_impl::parens( const int & noden,
					    std::ostream & o) const
  /*!
    Outputs the marginal tree in Newick format
  */
  {
    double time;

    if( left[noden] == -1 ) 
      {
	assert( (mi+((mi+noden)->abv))->time >= 0. );
	o << noden+1 << ':' <<(mi+((mi+noden)->abv))->time;
      }
    else
      {
	o << '(';
	parens( left[noden], o ) ;
	o << ',';
	parens( right[noden], o ) ;
	if( (mi+noden)->abv == -1 ) 
	  {
	    o << ");";
	  }
	else 
	  {
	    time = (mi + (mi+noden)->abv )->time - (mi+noden)->time ;
	    assert(time >= 0.);
	    o << "):"<<time;
	  }
      }
    return o;
  }

  void newick_stream_marginal_tree_impl::init()
  /*!
    Initializes the left and right vectors
  */
  {
    for(int i=0;i<2*nsam-2;++i)
      {
	if( left[ (mi+i)->abv ] == -1 ) 
	  {
	    left[(mi+i)->abv] = i;
	  }
	else 
	  {
	    right[ (mi+i)->abv] = i;
	  }
      }

  }

  newick_stream_marginal_tree::newick_stream_marginal_tree( const marginal & m ) :
    impl( new newick_stream_marginal_tree_impl(m) )
  {
  }

  newick_stream_marginal_tree::newick_stream_marginal_tree( const marginal * m ) :
    impl( new newick_stream_marginal_tree_impl(m) )
  {
  }

  newick_stream_marginal_tree::newick_stream_marginal_tree( arg::const_iterator m ) :
    impl( new newick_stream_marginal_tree_impl(m) )
  {
  }

  newick_stream_marginal_tree::newick_stream_marginal_tree( arg::iterator m ) :
    impl( new newick_stream_marginal_tree_impl(m) )
  {
  }

  newick_stream_marginal_tree::~newick_stream_marginal_tree( )
  {
    delete impl;
  }

  std::vector<node> newick_stream_marginal_tree::get_tree() const
  /*!
    if a tree has been read in from a stream, return it, else
    return an empty tree
  */
  {
    return impl->tree;
  }

  std::ostream & newick_stream_marginal_tree::print( std::ostream & o ) const
  /*!
    Write the marginal tree in Newick format to stream o
  */
  {
    impl->parens( 2*(impl->nsam)-2,o);
    return o;
  }

  std::istream & newick_stream_marginal_tree::read( std::istream & i )
  /*!
    Read a Newick tree from stream i
    \warning Not implemented!
  */
  {
    return i;
  }

  std::ostream & operator<<(std::ostream & o, const newick_stream_marginal_tree & n)
  /*!
    \return n.print(o);
  */
  {
    return n.print(o);
  }

  std::istream & operator>>(std::istream & i, newick_stream_marginal_tree & n)
  /*!
    \return n.read(o);
  */
  {
    return n.read(i);
  }
}
