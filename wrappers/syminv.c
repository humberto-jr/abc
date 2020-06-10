#include "wrappers.h"

void syminv_(double a[], const int lda, const int n, int *ierr)
{
	init_gpu();

	/* NOTE: both upper and lower parts are used. */
	for (int i = 0; i < n; ++i)
	{
		for (int j = (i + 1); j < n; ++j)
		{
			a[j*n + i] = a[i*n + j];
		}
	}

	magma_int_t *ipiv = NULL;
	magma_imalloc_cpu(&ipiv, n);

	double *a_gpu = NULL;
	magma_dmalloc(&a_gpu, n*n);

	magma_setmatrix(n, n, sizeof(double),
	                a, lda, a_gpu, lda, gpu_queue);

	magma_int_t info = 1;
	magma_dgetrf_gpu(n, n, a_gpu, lda, ipiv, &info);

	if (info != 0)
	{
		PRINT_ERROR("magma_dgetrf_gpu() failed with error code %d\n", info)
		*ierr = (int) info;
		return;
	}

	const magma_int_t lwork = n*magma_get_dgetri_nb(n);

	double *work = NULL;
	magma_dmalloc(&work, lwork);

	info = 1;
	magma_dgetri_gpu(n, a_gpu, n, ipiv, work, lwork, &info);

	if (info != 0)
	{
		PRINT_ERROR("magma_dgetri_gpu() failed with error code %d\n", info)
		*ierr = (int) info;
		return;
	}

	magma_getmatrix(n, n, sizeof(double),
	                a_gpu, lda, a, lda, gpu_queue);

	magma_free(ipiv);
	magma_free(work);
	magma_free(a_gpu);

	*ierr = (int) info;
}
