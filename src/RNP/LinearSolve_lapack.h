#ifndef _RNP_LINEAR_SOLVE_LAPACK_H_
#define _RNP_LINEAR_SOLVE_LAPACK_H_

#include "LinearSolve.h"
#include <cstdlib>

//#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifndef RNP_FORTRAN_NAME
# define RNP_FORTRAN_NAME(LCASE,UCASE) F77_FUNC(LCASE,UCASE)
#endif

extern "C" void RNP_FORTRAN_NAME(zgesv,ZGESV)(const long int &n, const long int &nrhs, std::complex<double> *a, 
	const long int &lda, long int *ipiv, std::complex<double> *b, const long int &ldb, long int *info);

extern "C" void RNP_FORTRAN_NAME(zgesvx,ZGESVX)(const char *fact, const char *trans,
	const long int &n, const long int &nrhs,
	std::complex<double> *a, const long int &lda,
	std::complex<double> *af, const long int &ldaf,
	long int *ipiv,
	char *equed, double *r, double *c,
	std::complex<double> *b, const long int &ldb,
	std::complex<double> *x, const long int &ldx,
	double *rcond, double *ferr, double *berr,
	std::complex<double> *work, double *rwork, long int *info);

namespace RNP{

template <>
template <>
inline LinearSolve<'N'>::LinearSolve(size_t n, size_t nRHS, std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, int *info, size_t *pivots){
	if(NULL != info){ *info = 0; }
	if(0 == n || nRHS == 0){ return; }
	
	size_t *ipiv = pivots;
	if(NULL == pivots){
		ipiv = (size_t*)malloc(sizeof(size_t)*n);
	}
	
	long int ret;
	zgesv_(n, nRHS, a, lda, (long int*)ipiv, b, ldb, &ret);
	
	if(NULL == pivots){
		free(ipiv);
	}
	if(NULL != info){ *info = ret; }
}
/*
inline LinearSolvePrecise(size_t n, size_t nRHS, std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, int *info, size_t *pivots){
	if(NULL != info){ *info = 0; }
	if(0 == n || nRHS == 0){ return; }
	
	size_t *ipiv = pivots;
	if(NULL == pivots){
		ipiv = (size_t*)malloc(sizeof(size_t)*n);
	}
	
	std::complex<double> *af = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n*n);
	double *rwork = (double*)malloc(sizeof(double) * 4*n);
	double *r = rwork + 2*n;
	double *c = r + n;
	
	char equed = "B";
	long int ret;
	zgesvx_(
		"E", "N", n, nrhs,
		a, lda, af, n,
		ipiv,
		equed, r, c,
		b, ldb,
		std::complex<double> *x, const long int &ldx,
		double *rcond, double *ferr, double *berr,
		std::complex<double> *work, double *rwork, long int *info);
	zgesv_(n, nRHS, a, lda, (long int*)ipiv, b, ldb, &ret);
	
	free(af);
	
	if(NULL == pivots){
		free(ipiv);
	}
	if(NULL != info){ *info = ret; }
}
*/

}; // namespace RNP

#endif // _RNP_LINEAR_SOLVE_LAPACK_H_
