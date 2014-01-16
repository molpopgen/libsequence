#ifndef __SEQ_SAMFUNCTIONS_HPP__
#define __SEQ_SAMFUNCTIONS_HPP__

#include <Sequence/samrecord.hpp>

namespace Sequence
{
  unsigned alignment_length( const samrecord & b );
  unsigned insertion_distance( const samrecord & b );
  unsigned deletion_distance( const samrecord & b );
  unsigned ngaps( const samrecord & b );
  unsigned mismatches( const samrecord & b );
}
#endif
