#include "wrappers.h"

void DSYR_(const char uplo,
           const int n,
           const double alpha,
           const double x[],
           const int incx,
           double a[],
           const int lda)
{
	init_gpu();

	const magma_uplo_t a_form = (uplo == 'u'? MagmaUpper : MagmaLower);

	double *x_gpu = NULL, *a_gpu = NULL;

	magma_dmalloc(&x_gpu, n);
	magma_dmalloc(&a_gpu, n*n);

	magma_setvector_async(n, sizeof(double), x, incx, x_gpu, incx, gpu_queue);
	magma_setmatrix_async(n, n, sizeof(double), a, lda, a_gpu, lda, gpu_queue);

	magma_queue_sync(gpu_queue);
	magma_dsyr(a_form, n, alpha, x_gpu, incx, a_gpu, lda, gpu_queue);

	magma_getmatrix(n, n, sizeof(double), a_gpu, lda, a, lda, gpu_queue);

	magma_free(x_gpu);
	magma_free(a_gpu);
}
