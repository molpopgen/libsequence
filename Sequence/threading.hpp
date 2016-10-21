#ifndef SEQUENCE_THREADING_HPP_
#define SEQUENCE_THREADING_HPP_

#include <memory>
#include <tbb/task_scheduler_init.h>

/*!
 * \defgroup threads Parallel computing
 * \brief Classes & functions related to using parallelized algorithms.
 *
 * libsequence uses the Intel TBB library for parallel processing.
 * With TBB, the number of threads is determined automatically. The
 * default/automatic behavior can be changed by creating an instance
 * of a tbb::task_scheduler_init object:
 *
 * \code
 * #include <tbb/task_scheduler_init.h>
 *
 * //restrict parallelism to max_threads
 * tbb::task_scheduler_init init(max_threads);
 * \endcode
 *
 * Doing the above requires that your program further link to TBB,
 * via -ltbb. When using the default scheme (e.g., your program does
 * not manage its own task_scheduler_init), you do not need to link to 
 * TBB, as libsequence itself is linked.
 *
 * For cases where libsequence is used as a back-end for projects in
 * other languages like R or Python, we provide Sequence::init_tbb which
 * returns a std::unique_ptr<tbb::task_scheduler_init>.  This function can be
 * wrapped, allowing control over threading as desired.
 *
 * Unfortunately, Doxygen limits functions to belonging to a single group,
 * making it difficult to auto-generate a list of parallelized functions.
 * The following list of threaded functions in other groups is manually
 * maintained, and may thus be incomplete:
 *
 * - Sequence::nSL_t
 * - Sequence::snSL
 * - Sequence::GarudStats
 * - Sequence::Recombination::Disequilibrium
 *
 * While there are lots of calculations that may be parallelized, not all
 * will benefit.  We have to test each function individually and make the 
 * determination.  Fortunately, TBB allows us to do this in the implementation
 * without changing the public interface, meaning we can take our time and slowly
 * introduce more parallelism into the library.
 */

namespace Sequence
{
	/*! \brief Return a unique_ptr to a tbb::task_scheduler_init
	 * 
	 * This exists mainly for libsequence-based extensions in other languages.
	 * This function is not necessary in a standalone program where you can rely
	 * on creating a tbb::task_scheduler_init on the stack.
	 *
	 * \return std::unique_ptr<tbb::task_scheduler_init>
	 *
	 * \ingroup threads
	 */
    std::unique_ptr<tbb::task_scheduler_init>
    init_tbb(int nthreads = tbb::task_scheduler_init::automatic);
}

#endif
