#include "wrappers.h"

void DSYTRF_(const char uplo,
             const int n,
             double a[],
             const int lda,
             int ipiv[],
             double work[],
             const int lwork,
             int *info)
{
	init_gpu();

	const magma_uplo_t fill_mode = (uplo == 'u'? MagmaUpper : MagmaLower);

	double *a_gpu = NULL;

	magma_dmalloc(&a_gpu, lda*n);

	magma_setmatrix(n, n, sizeof(double), a, lda, a_gpu, lda, gpu_queue);

	magma_dsytrf(fill_mode, n, a_gpu, lda, ipiv, info);

	if (*info != 0)
	{
		PRINT_ERROR("magma_dsytrf() failed with error code %d\n", *info)
		exit(EXIT_FAILURE);
	}

	magma_getmatrix(n, n, sizeof(double), a_gpu, lda, a, lda, gpu_queue);

	magma_free(a_gpu);
}
