#if !defined(CLIB_HEADER)
	#define CLIB_HEADER

	/* C89 to C95 */
	#include <math.h>
	#include <time.h>
	#include <float.h>
	#include <wchar.h>
	#include <errno.h>
	#include <stdio.h>
	#include <ctype.h>
	#include <assert.h>
	#include <iso646.h>
	#include <wctype.h>
	#include <limits.h>
	#include <locale.h>
	#include <setjmp.h>
	#include <signal.h>
	#include <stdarg.h>
	#include <stddef.h>
	#include <stdlib.h>
	#include <string.h>

	/* C POSIX library */
	#include <unistd.h>

	/* C99 */
	#if (__STDC_VERSION__ >= 199901L)
		#include <fenv.h>
		#include <iso646.h>
		#include <stdint.h>
		#include <tgmath.h>
		#include <complex.h>
		#include <stdbool.h>
		#include <inttypes.h>
	#endif

	/* C11 and beyond */
	#if (__STDC_VERSION__ > 199901L)
		#include <uchar.h>
		#include <threads.h>
		#include <stdalign.h>
		#include <stdatomic.h>
		#include <stdnoreturn.h>
	#endif

	/* NOTE: just a few types, to complete as needed. */
	enum c_type
	{
		type_int,
		type_char,
		type_bool,
		type_float,
		type_double,
		type_long_int,
		type_long_double,
		type_unsigned_int,
		type_unsigned_char,
		type_unknown
	};

	typedef enum c_type c_type;
#endif
