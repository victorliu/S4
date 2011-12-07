#ifndef _RNP_EIGENSYSTEMS_LAPACK_H_
#define _RNP_EIGENSYSTEMS_LAPACK_H_

#include <cstddef>
#include <complex>

//#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifndef RNP_FORTRAN_NAME
# define RNP_FORTRAN_NAME(LCASE,UCASE) F77_FUNC(LCASE,UCASE)
#endif

#include <essl.h>

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
	int iopt = 0;
	if(vl != NULL){ iopt = 1; }
	if(vr != NULL){ iopt = 1; }

	if((size_t)-1 == lwork_){
		work_[0] = 0;
		return 0;
	}
	 int    esvzgeev   (int, void *, int, _ESVCOM *, void *, int, const int *, int, double *, int);

	zgeev(iopt, a, lda, eval, vr, ldvr, NULL, n, rwork_, 0);
	
	return 0;
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
}; // namespace RNP

#endif // _RNP_EIGENSYSTEMS_LAPACK_H_
