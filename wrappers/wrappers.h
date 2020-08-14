#if !defined(WRAPPERS_HEADER)
	#define WRAPPERS_HEADER

	#include "c_lib.h"
	#include <magma_v2.h>
	#include <magma_lapack.h>

	static magma_queue_t gpu_queue;

	/******************************************************************************

	 Macro PRINT_ERROR(): writes an error message in the C stderr using the format
	 "# file.c, function(), line n: formatted message here".

	******************************************************************************/

	#define PRINT_ERROR(format, ...)                                          \
	{                                                                         \
	  fprintf(stderr, "# %s, %s(), line %d: ", __FILE__, __func__, __LINE__); \
	  fprintf(stderr, format, ##__VA_ARGS__);                                 \
	}

	/******************************************************************************

	 Function max(): return the max between two integers a and b.

	******************************************************************************/

	inline static int max(const int a, const int b)
	{
		return (a > b? a : b);
	}

	/******************************************************************************

	 Function init_gpu(): allocate resources needed by MAGMA. Invoked only once.

	******************************************************************************/

	static inline void init_gpu()
	{
		static bool not_ready = true;

		if (not_ready)
		{
			if (magma_init() != MAGMA_SUCCESS)
			{
				PRINT_ERROR("%s\n", "magma_init() failed")
				exit(EXIT_FAILURE);
			}

			magma_queue_create(0, &gpu_queue);
			not_ready = false;
		}
	}
#endif
