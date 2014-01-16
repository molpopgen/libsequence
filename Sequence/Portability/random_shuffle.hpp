#ifndef __SEQUENCE_PORTABILITY_RANDOM_SHUFFLE_HPP__
#define __SEQUENCE_PORTABILITY_RANDOM_SHUFFLE_HPP__
#include <iterator>
#include <boost/ref.hpp>

namespace Sequence
{
  template<typename iterator,typename generator>
  inline void random_shuffle( iterator __first, iterator __last, generator & __g )
  {
    typedef typename std::iterator_traits<iterator>::difference_type __dt;
    if(__first==__last)return;
    for( iterator __i = __first + 1 ; __i != __last ; ++__i )
      {  
	std::iter_swap(__i,__first+__dt(__g(boost::cref(((__i -__first)+1)))));
      }
  }

  template<typename iterator, typename generator>
  inline void random_shuffle(iterator __first, iterator __last, const generator & __g)
  {
    typedef typename std::iterator_traits<iterator>::difference_type __dt;
    if(__first==__last)return;
    for( iterator __i = __first + 1 ; __i != __last ; ++__i )
      {  
	std::iter_swap(__i,__first+__dt(__g(boost::cref(((__i -__first)+1)))));
      }
  }
}

#endif
