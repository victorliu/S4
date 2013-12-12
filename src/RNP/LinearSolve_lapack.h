#ifndef _RNP_LINEAR_SOLVE_LAPACK_H_
#define _RNP_LINEAR_SOLVE_LAPACK_H_

#include "LinearSolve.h"
#include <cstdlib>
#include <stdint.h>

//#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifndef RNP_FORTRAN_NAME
# define RNP_FORTRAN_NAME(LCASE,UCASE) F77_FUNC(LCASE,UCASE)
#endif

typedef int integer;

extern "C" void RNP_FORTRAN_NAME(zgesv,ZGESV)(const integer &n, const integer &nrhs, std::complex<double> *a, 
	const integer &lda, integer *ipiv, std::complex<double> *b, const integer &ldb, integer *info);

extern "C" void RNP_FORTRAN_NAME(zgesvx,ZGESVX)(const char *fact, const char *trans,
	const integer &n, const integer &nrhs,
	std::complex<double> *a, const integer &lda,
	std::complex<double> *af, const integer &ldaf,
	integer *ipiv,
	char *equed, double *r, double *c,
	std::complex<double> *b, const integer &ldb,
	std::complex<double> *x, const integer &ldx,
	double *rcond, double *ferr, double *berr,
	std::complex<double> *work, double *rwork, integer *info);

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
	
	integer ret;
	RNP_FORTRAN_NAME(zgesv,ZGESV)(n, nRHS, a, lda, (integer*)ipiv, b, ldb, &ret);
	
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
	integer ret;
	zgesvx_(
		"E", "N", n, nrhs,
		a, lda, af, n,
		ipiv,
		equed, r, c,
		b, ldb,
		std::complex<double> *x, const integer &ldx,
		double *rcond, double *ferr, double *berr,
		std::complex<double> *work, double *rwork, integer *info);
	zgesv_(n, nRHS, a, lda, (integer*)ipiv, b, ldb, &ret);
	
	free(af);
	
	if(NULL == pivots){
		free(ipiv);
	}
	if(NULL != info){ *info = ret; }
}
*/

}; // namespace RNP

#endif // _RNP_LINEAR_SOLVE_LAPACK_H_
