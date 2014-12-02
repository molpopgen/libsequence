#ifndef __SEQ_SAMFUNCTIONS_HPP__
#define __SEQ_SAMFUNCTIONS_HPP__

#include <Sequence/samrecord.hpp>
#include <Sequence/bamrecord.hpp>

namespace Sequence
{
  unsigned alignment_length( const samrecord & b );
  unsigned insertion_distance( const samrecord & b );
  unsigned deletion_distance( const samrecord & b );
  unsigned ngaps( const samrecord & b );
  unsigned mismatches( const samrecord & b );

#ifdef HAVE_HTSLIB
  unsigned alignment_length( const bamrecord & b );
  unsigned insertion_distance( const bamrecord & b );
  unsigned deletion_distance( const bamrecord & b );
  unsigned ngaps( const bamrecord & b );
  unsigned mismatches( const bamrecord & b );
#endif
}
#endif
