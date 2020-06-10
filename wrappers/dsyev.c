#include "wrappers.h"

void dsyev_(const char jobz,
            const char uplo,
            const int n,
            double a[],
            const int lda,
            double w[])
{
	init_gpu();

	const magma_vec_t job_mode
		= (jobz == 'n'? MagmaNoVec : MagmaVec);

	const magma_uplo_t fill_mode
		= (uplo == 'u'? MagmaUpper : MagmaLower);

	const magma_int_t nb
		= magma_get_dsytrd_nb((magma_int_t) n);

	const magma_int_t lwork
		= (jobz == 'n'? 2*n + n*nb : max(2*n + n*nb, 1 + 6*n + 2*n*n));

	const magma_int_t liwork
		= (jobz == 'n'? n : 3 + 5*n);

	double *a_gpu = NULL, *work = NULL;

	magma_dmalloc(&a_gpu, n*n);
	magma_dmalloc_cpu(&work, lwork);

	magma_int_t *iwork = NULL, info = 1;

	magma_imalloc_cpu(&iwork, liwork);

	magma_dsetmatrix(n, n, a, lda, a_gpu, lda, gpu_queue);

	magma_dsyevd_gpu(job_mode, fill_mode, n, a_gpu, lda,
	                 w, a, n, work, lwork, iwork, liwork, &info);

	if (info != 0)
	{
		PRINT_ERROR("magma_dsyevd_gpu() failed with error code %d\n", info)
		exit(EXIT_FAILURE);
	}

	if (jobz != 'n')
	{
		magma_dgetmatrix(n, n, a_gpu, lda, a, lda, gpu_queue);
	}

	magma_free(a_gpu);

	magma_free_cpu(work);
	magma_free_cpu(iwork);
}
