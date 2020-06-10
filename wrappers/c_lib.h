#if !defined(CLIB_HEADER)
	#define CLIB_HEADER

	/* C89 to C99 */
	#include <math.h>
	#include <time.h>
	#include <fenv.h>
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
	#include <stdint.h>
	#include <stdlib.h>
	#include <string.h>
	#include <tgmath.h>
	#include <stdbool.h>
	#include <complex.h>
	#include <inttypes.h>

	/* C POSIX library */ 
	#include <unistd.h>

	#if (__STDC_VERSION__ > 199901L)
		#include <uchar.h>
		#include <threads.h>
		#include <stdalign.h>
		#include <stdatomic.h>
		#include <stdnoreturn.h>
	#endif
#endif
