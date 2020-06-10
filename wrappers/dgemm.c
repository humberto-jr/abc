#include "wrappers.h"

void dgemm_(const char trans_a,
            const char trans_b,
            const int m,
            const int n,
            const int k,
            const double alpha,
            const double a[],
            const int lda,
            const double b[],
            const int ldb,
            const double beta,
            double c[],
            const int ldc)
{
	init_gpu();

	const magma_trans_t a_form
		= (trans_a == 'n'? MagmaNoTrans : MagmaTrans);

	const magma_trans_t b_form
		= (trans_b == 'n'? MagmaNoTrans : MagmaTrans);

	double *a_gpu = NULL, *b_gpu = NULL, *c_gpu = NULL;

	magma_dmalloc(&a_gpu, m*k);
	magma_dmalloc(&b_gpu, k*n);
	magma_dmalloc(&c_gpu, m*n);

	magma_setmatrix_async(m, k, sizeof(double), a, lda, a_gpu, lda, gpu_queue);
	magma_setmatrix_async(k, n, sizeof(double), b, ldb, b_gpu, ldb, gpu_queue);

	if (beta != 0.0)
	{
		magma_setmatrix_async(m, n, sizeof(double), c, ldc, c_gpu, ldc, gpu_queue);
	}

	magma_queue_sync(gpu_queue);

	magmablas_dgemm(a_form, b_form, m, n, k,
	                alpha, a_gpu, lda, b_gpu, ldb, beta, c_gpu, ldc, gpu_queue);

	magma_getmatrix(m, n, sizeof(double), c_gpu, ldc, c, ldc, gpu_queue);

	magma_free(a_gpu);
	magma_free(b_gpu);
	magma_free(c_gpu);
}
