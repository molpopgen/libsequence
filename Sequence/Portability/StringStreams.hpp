/*! \file StringStreams.hpp Portability header
  This header is part of the libsequence package.

  Typedefs are provided in namespace Sequence for the 
  standard string stream headers.  By default, the 
  ANSI/ISO sstream header is preferred, and is included
  either if no macros are defined or if HAVE_SSTREAM is 
  defined.  

  To use the "traditional" header strstream, 
  define HAVE_STRSTREAM and leave HAVE_SSTREAM undefined.
*/

#ifndef __STRING_STREAMS_HPP__
#define __STRING_STREAMS_HPP__

#if defined(HAVE_STRSTREAM) & (!(defined(HAVE_SSTREAM)))
#include <strstream>
#else
#include <sstream>
#endif

namespace Sequence
{
#if defined(HAVE_STRSTREAM) && (!(defined(HAVE_SSTREAM)))
#warning "using deprecated (char * based) string streams"
  typedef std::istrstream istr;
  typedef std::ostrstream ostr;
#else
  typedef std::istringstream istr;
  typedef std::ostringstream ostr;
#endif
}

#endif
