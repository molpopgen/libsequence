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

#ifndef __ENSURE_FLOATING_HPP__
#define __ENSURE_FLOATING_HPP__

#include <Sequence/preferFloatingTypes.hpp>
#include <boost/mpl/identity.hpp>

namespace Sequence
{
  template <typename T, typename T2, typename floatingType = double>
  struct ensureFloating
  /*! \struct ensureFloating Sequence/ensureFloating.hpp
    \ingroup metaprogramming
    A metaprogram to ensure that a floating type is chosen.  During 
    compilation, if either T or T2 are floating types, the floating
    type is selected according to the rules of Sequence::preferFloatingTypes.
    If neither T nor T2 is a floating type, floatingType is selected.
    Compile-time assertions check that T and T2 are arithmetic types,
    and the floatingType is indeed a floating type.
   */
  {
    BOOST_STATIC_ASSERT( (boost::is_arithmetic<T>::value) &&
			 (boost::is_arithmetic<T2>::value) );
    BOOST_STATIC_ASSERT( (boost::is_float<floatingType>::value) );
    /*!
      resolves to the selected type
    */
    typedef typename boost::mpl::if_<
      typename boost::mpl::not_<
      typename boost::is_float<typename preferFloatingTypes<T,T2>::type>
    >,
      typename boost::mpl::identity<floatingType>::type,
      typename preferFloatingTypes<T,T2>::type >::type type;
  };
}
#endif
