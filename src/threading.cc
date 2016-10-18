#include <Sequence/threading.hpp>

namespace Sequence
{
    std::unique_ptr<tbb::task_scheduler_init>
    init_tbb(int nthreads)
	{
		return std::unique_ptr<tbb::task_scheduler_init>(new tbb::task_scheduler_init(nthreads));
	}
}
