#ifndef _RNP_EIGENSYSTEMS_LAPACK_H_
#define _RNP_EIGENSYSTEMS_LAPACK_H_

#include <cstddef>
#include <complex>

//#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifndef RNP_FORTRAN_NAME
# define RNP_FORTRAN_NAME(LCASE,UCASE) F77_FUNC(LCASE,UCASE)
#endif

typedef int integer;

extern "C" void RNP_FORTRAN_NAME(zgeev,ZGEEV)(const char *jobvl, const char *jobvr, const integer &n, 
	std::complex<double> *a, const integer &lda, std::complex<double> *w, std::complex<double> *vl, 
	const integer &ldvl, std::complex<double> *vr, const integer &ldvr, std::complex<double> *work, 
	const integer &lwork, double *rwork, integer *info);

extern "C" void RNP_FORTRAN_NAME(zgees,ZGEES)(const char *jobvs, const char *sort, integer (*select)(), const integer &n,
	std::complex<double> *a, const integer &lda, const integer &sdim, std::complex<double> *w,
	std::complex<double> *vs, const integer &ldvs, std::complex<double> *work, const integer &lwork,
	double *rwork, integer *bwork, integer *info);

extern "C" void RNP_FORTRAN_NAME(dgeev,DGEEV)(const char *jobvl, const char *jobvr, const integer &n, 
	double *a, const integer &lda, double *wr, double *wi, double *vl, 
	const integer &ldvl, double *vr, const integer &ldvr, double *work, 
	const integer &lwork, integer *info);

#include <iostream>

namespace RNP{

// Eigensystem computes for an N-by-N complex nonsymmetric matrix A, the
// eigenvalues and, optionally, the left and/or right eigenvectors.
//
// The right eigenvector v(j) of A satisfies
//                  A * v(j) = lambda(j) * v(j)
// where lambda(j) is its eigenvalue.
// The left eigenvector u(j) of A satisfies
//               u(j)**H * A = lambda(j) * u(j)**H
// where u(j)**H denotes the conjugate transpose of u(j).
//
// The computed eigenvectors are normalized to have Euclidean norm
// equal to 1 and largest component real.
//
// Arguments
// =========
//
// n       The order of the matrix A. n >= 0.
//
// a       (input/output) dimension (lda,n)
//         On entry, the N-by-N matrix A.
//         On exit, A has been overwritten.
//
// lda     The leading dimension of the array A.  lda >= n.
//
// eval    (output) dimension (n)
//         eval contains the computed eigenvalues.
//
// vl      (output) dimension (ldvl,n)
//         If not NULL, the left eigenvectors u(j) are stored one
//         after another in the columns of vl, in the same order
//         as their eigenvalues.
//         u(j) = vl(:,j), the j-th column of vl.
//
// ldvl    The leading dimension of the array vl. If vl != NULL, ldvl >= N.
//
// vr      (output) dimension (ldvr,n)
//         If not NULL, the right eigenvectors v(j) are stored one
//         after another in the columns of vr, in the same order
//         as their eigenvalues.
//         v(j) = vr(:,j), the j-th column of vr.
//
// ldvr    The leading dimension of the array vr. vr != NULL, ldvr >= n.
//
// work    dimension 2*N
// rwork   dimension 2*N
//
// return  = 0:  successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value.
//         > 0:  if INFO = i, the QR algorithm failed to compute all the
//               eigenvalues, and no eigenvectors have been computed;
//               elements and i+1:N of W contain eigenvalues which have
//               converged.
inline int Eigensystem(size_t n,
	std::complex<double> *a, size_t lda,
	std::complex<double> *eval,
	std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
	std::complex<double> *work_, double *rwork_, size_t lwork_)
{
	if(n == 0) {
		return 0;
	}
	char jobvl[2] = "N";
	char jobvr[2] = "N";
	if(vl != NULL){ jobvl[0] = 'V'; }
	if(vr != NULL){ jobvr[0] = 'V'; }
	integer info;

	if((size_t)-1 == lwork_){
		RNP_FORTRAN_NAME(zgeev,ZGEEV)(jobvl, jobvr, n, a, lda, eval, vl, ldvl, vr, ldvr, work_, -1, rwork_, &info);
		return 0;
	}
	
	std::complex<double> *work = work_;
	double *rwork = rwork_;
	if(NULL == rwork_){
		rwork = new double[2*n];
	}
	if(0 == lwork_){ lwork_ = 2*n; }
	integer lwork = lwork_;
	if(NULL == work_ || lwork < (integer)(2*n)){
		lwork = -1;
		std::complex<double> zlen;
		RNP_FORTRAN_NAME(zgeev,ZGEEV)(jobvl, jobvr, n, a, lda, eval, vl, ldvl, vr, ldvr, &zlen, lwork, rwork, &info);
		lwork = (integer)(zlen.real());
		work = new std::complex<double>[lwork];
	}
	RNP_FORTRAN_NAME(zgeev,ZGEEV)(jobvl, jobvr, n, a, lda, eval, vl, ldvl, vr, ldvr, work, lwork, rwork, &info);
	
	if(NULL == work_ || lwork < (integer)(2*n)){
		delete [] work;
	}
	if(NULL == rwork_){
		delete [] rwork;
	}
	
	return info;
}


inline int Eigensystem_real(size_t n,
	double *a, size_t lda,
	std::complex<double> *eval,
	std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
	double *work_, size_t lwork_)
{
	if(n == 0) {
		return 0;
	}
	char jobvl[2] = "N";
	char jobvr[2] = "N";
	if(vl != NULL){ jobvl[0] = 'V'; }
	if(vr != NULL){ jobvr[0] = 'V'; }
	integer info;
	
	double *peval = (double*)eval;

	if((size_t)-1 == lwork_){
		RNP_FORTRAN_NAME(dgeev,DGEEV)(jobvl, jobvr, n, a, lda, peval, peval+n, (double*)vl, 2*ldvl, (double*)vr, 2*ldvr, work_, -1, &info);
		return 0;
	}
	
	double *work = work_;
	if(0 == lwork_){ lwork_ = 4*n; }
	integer lwork = lwork_;
	if(NULL == work_ || lwork < (integer)(4*n)){
		lwork = -1;
		double zlen;
		RNP_FORTRAN_NAME(dgeev,DGEEV)(jobvl, jobvr, n, a, lda, peval, peval+n, (double*)vl, 2*ldvl, (double*)vr, 2*ldvr, &zlen, lwork, &info);
		lwork = (integer)zlen;
		work = new double[lwork];
	}
	RNP_FORTRAN_NAME(dgeev,DGEEV)(jobvl, jobvr, n, a, lda, peval, peval+n, (double*)vl, 2*ldvl, (double*)vr, 2*ldvr, work, lwork, &info);
	
	double *wi = work+2*n;
	memcpy(wi, peval+n, sizeof(double)*n);
	for(size_t j = 0; j < n; ++j){
		memcpy(&a[j*lda], &vr[j*ldvr], sizeof(double)*n);
	}
	size_t j = n;
	while(j --> 0){
		eval[j] = std::complex<double>(peval[j], wi[j]);
		if(0 == wi[j]){
			if(NULL != vl){
				size_t i = n;
				while(i --> 0){
					vl[i+j*ldvl] = *(((double*)vl)+i+j*2*ldvl);
				}
			}
			if(NULL != vr){
				size_t i = n;
				while(i --> 0){
					vr[i+j*ldvr] = *(((double*)vr)+i+j*2*ldvr);
				}
			}
		}else{
			if(j == 0){
				break;
			}
			eval[j-1] = std::complex<double>(peval[j-1], wi[j-1]);
			if(NULL != vl){
				// TODO
			}
			if(NULL != vr){
				for(size_t i = 0; i < n; ++i){
					vr[i+(j-1)*ldvr] = std::complex<double>(a[i+(j-1)*lda],  a[i+j*lda]);
					vr[i+ j   *ldvr] = std::complex<double>(a[i+(j-1)*lda], -a[i+j*lda]);
				}
			}
			j--;
		}
	}
	
	if(NULL == work_ || lwork < (integer)(4*n)){
		delete [] work;
	}
	
	return info;
}

/*
// zheev
template <char uplo>
int HermitianEigensystem(size_t n, 
	std::complex<double> *a, size_t lda,
	std::complex<double> *eval,
	std::complex<double> *work, double *rwork);

// dsyev
template <char uplo>
int SymmetricEigensystem(size_t n, 
	double *a, size_t lda,
	double *eval,
	double *work);
*/
} // namespace RNP

#endif // _RNP_EIGENSYSTEMS_LAPACK_H_
