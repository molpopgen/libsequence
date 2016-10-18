#ifndef SEQUENCE_THREADING_HPP_
#define SEQUENCE_THREADING_HPP_

#include <memory>
#include <tbb/task_scheduler_init.h>

/*!
 * \defgroup threads Parallel computing
 * \brief Classes & functions related to using parallelized algorithms.
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
	 */
    std::unique_ptr<tbb::task_scheduler_init>
    init_tbb(int nthreads = tbb::task_scheduler_init::automatic);
}

#endif
