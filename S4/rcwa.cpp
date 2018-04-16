/* Copyright (C) 2009-2011, Stanford University
 * This file is part of S4
 * Written by Victor Liu (vkl@stanford.edu)
 *
 * S4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * S4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "config.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <float.h>
#include "rcwa.h"
#include "fmm/fft_iface.h"
#include <TBLAS.h>
#ifdef HAVE_BLAS
# include <TBLAS_ext.h>
#endif
#include <LinearSolve.h>
#ifdef HAVE_LAPACK
# include <LinearSolve_lapack.h>
# include <Eigensystems_lapack.h>
#else
# include <Eigensystems.h>
#endif

#ifdef DUMP_MATRICES
# define RNP_OUTPUT_MATHEMATICA
# define DUMP_STREAM std::cerr
//# define DUMP_STREAM (omega.real() > 1.91637 ? std::cerr : std::cout)
# include <IO.h>
#endif
#include <cstdio>

#include <numalloc.h>

static inline void* rcwa_malloc(size_t size){
	void *ret = malloc_aligned(size, 16);
	//memset(ret, 0x0, size);
	return ret;
}
static inline void rcwa_free(void *ptr){
	free_aligned(ptr);
}

#ifdef HAVE_LAPACK
	extern "C" void zgelss_(
		const integer &m, const integer &n, const integer &nRHS,
		std::complex<double> *a, const integer &lda,
		std::complex<double> *b, const integer &ldb,
		double *s, const double &rcond, integer *rank,
		std::complex<double> *work, const integer &lwork,
		double *rwork, integer *info
	);
#endif

static void SingularLinearSolve(
	size_t m, size_t n, size_t nRHS, std::complex<double> *a, size_t lda,
	std::complex<double> *b, size_t ldb, const double &rcond
){
#ifdef HAVE_LAPACK
	integer info, rank;
	std::complex<double> dummy;
	zgelss_(m, n, nRHS, a, lda, b, ldb, NULL, rcond, &rank, &dummy, -1, NULL, &info);
	integer lwork = (integer)dummy.real();
	std::complex<double> *work = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>) * lwork);
	double *rwork = (double*)rcwa_malloc(sizeof(double) * 6*m);

	zgelss_(m, n, nRHS, a, lda, b, ldb, rwork, rcond, &rank, work, lwork, rwork+m, &info);

	rcwa_free(rwork);
	rcwa_free(work);
#else
	RNP::LinearSolve<'N'>(n, nRHS, a, lda, b, ldb);
#endif
}


template <typename T>
static void Swap(int n, T* x, int incx, T* y, int incy){
	while(n --> 0){
		std::swap(*x, *y);
		x += incx; y += incy;
	}
}

template <typename A, typename T>
static void Axpy(int n, const A &alpha, const T *x, int incx, T *y, int incy){
	while(n --> 0){
		*y += alpha*(*x);
		x += incx;
		y += incy;
	}
}
template <typename A, typename T>
static void Axpy(int m, int n, const A &alpha, const T *a, int lda, T *b, int ldb){
	for(int j = 0; j < n; ++j){
		Axpy(m, alpha, a, 1, b, 1);
		a += lda;
		b += ldb;
	}
}
template <typename S, typename T>
static void Scale(int n, const S &s, T *v, int incv = 1){
	while(n --> 0){
		*v *= s;
		v += incv;
	}
}

template <typename T>
static void Copy(int n, const T *src, int incsrc, T *dst, int incdst){
	if(1 == incsrc && 1 == incdst){
		memcpy(dst, src, sizeof(T)*n);
	}else{
		for(int j = 0; j < n; ++j){
			*dst = *src;
			src += incsrc;
			dst += incdst;
		}
	}
}
template <typename T>
static void Copy(int m, int n, const T *src, int ldsrc, T *dst, int lddst){
	if(m == ldsrc && m == lddst){
		memcpy(dst, src, sizeof(T)*m*n);
	}else{
		for(int j = 0; j < n; ++j){
			memcpy(dst, src, sizeof(T)*m);
			src += ldsrc;
			dst += lddst;
		}
	}
}
template <typename A, typename T>
static void SetMatrix(int m, int n, const A &diag, const A &offdiag, T *a, int lda){
	for(int j = 0; j < n; ++j){
		int i;
		T *p = a;
		for(i = 0; i < j; ++i){
			*p++ = offdiag;
		}
		*p++ = diag;
		for(i++ ; i < m; ++i){
			*p++ = offdiag;
		}
		a += lda;
	}
}
template <typename T>
static void ZeroMatrix(int m, int n, T *a, int lda){
	if(m == lda){
		memset(a, 0, sizeof(T)*m*n);
	}else{
		for(int j = 0; j < n; ++j){
			memset(a, 0, sizeof(T)*m);
			a += lda;
		}
	}
}
extern "C" void zgemv_(
	const char *trans, const int &m, const int &n, const std::complex<double> &alpha,
	const std::complex<double> *a, const int &lda, const std::complex<double> *x, const int &incx,
	const std::complex<double> &beta, std::complex<double> *y, const int &incy
);
extern "C" void zgemm_(
	const char *transa, const char *transb, const int &m, const int &n, const int &k,
	const std::complex<double> &alpha, const std::complex<double> *a, const int &lda,
	const std::complex<double> *b, const int &ldb,
	const std::complex<double> &beta, std::complex<double> *c, const int &ldc
);
extern "C" void zgetrf_(
	const int &m, const int &n, std::complex<double> *a, const int &lda, int *ipiv, int *info
);
extern "C" void zgetrs_(
	const char *trans, const int &n, const int &nrhs, const std::complex<double> *a, const int &lda,
	const int *ipiv, std::complex<double> *b, const int &ldb, int *info
);
extern "C" void zgetri_(
	const int &n, const std::complex<double> *a, const int &lda,
	const int *ipiv, std::complex<double> *work, const int &lwork, int *info
);
static int Invert(size_t n, std::complex<double> *a, size_t lda, std::complex<double> *work, size_t lwork, size_t *iwork){
	int info;
	zgetrf_(n, n, a, lda, (int*)iwork, &info);
	zgetri_(n, a, lda, (int*)iwork, work, lwork, &info);
	return info;
}

static void Mult(size_t n, const double &alpha, const std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, const double &beta, std::complex<double> *c, size_t ldc){
	zgemm_("N","N", n, n, n, alpha, a, lda, b, ldb, beta, c, ldc);
}
static void Mult(size_t n, size_t m, const double &alpha, const std::complex<double> *a, size_t lda, std::complex<double> *b, const double &beta, std::complex<double> *c){
	zgemv_("N", n, m, alpha, a, lda, b, 1, beta, c, 1);
}
static int LUFactor(size_t n, std::complex<double> *a, size_t lda, size_t *ipiv){
	int info = 0;
	zgetrf_(n, n, a, lda, (int*)ipiv, &info); // A = P*L*U
	return info;
}
static int LUSolve(size_t n, size_t nRHS, const std::complex<double> *a, size_t lda, const size_t *ipiv, std::complex<double> *b, size_t ldb){
	int info;
	zgetrs_("N", n, nRHS, a, lda, (const int*)ipiv, b, ldb, &info);
	return info;
}
static void PrintMatrix(const char *name, size_t m, size_t n, const std::complex<double> *a, size_t lda){
	printf("%s = [\n", name);
	for(size_t i = 0; i < m; ++i){
		for(size_t j = 0; j < n; ++j){
			const std::complex<double> &val = a[i+j*lda];
			if(val.imag() < 0){
				printf(" % 18.16e-% 18.16ei", val.real(), -val.imag());
			}else{
				printf(" % 18.16e+% 18.16ei", val.real(), val.imag());
			}
		}
		printf("\n");
	}
	printf("];\n");
}

// kp = k-parallel matrix
// kp = omega^2 - Kappa = omega^2 - [  ky*epsinv*ky -ky*epsinv*kx ]
//                                  [ -kx*epsinv*ky  kx*epsinv*kx ]
//    = omega^2 + [ -ky ] [ epsinv ] [ ky -kx ]
//

static void MakeKPMatrix(
	std::complex<double> omega,
	size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv,
	int epstype,
	const std::complex<double> *kp_existing,
	std::complex<double> *kp,
	const size_t ldkp // leading dimension of kp, >= 2*n
){
	const size_t n2 = 2*n;
	if(NULL != kp_existing){
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, kp_existing,n2, kp,ldkp);
		return;
	}

	const std::complex<double> omega2 = omega*omega;
	if(EPSILON2_TYPE_BLKDIAG1_SCALAR == epstype || EPSILON2_TYPE_BLKDIAG2_SCALAR == epstype){
		RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., kp,ldkp);

		for(size_t i = 0; i < n; ++i){
			kp[(0+i)+(0+i)*ldkp] = omega2 - ky[i]*Epsilon_inv[0]*ky[i];
			kp[(n+i)+(0+i)*ldkp] = kx[i]*Epsilon_inv[0]*ky[i];
		}
		for(size_t i = 0; i < n; ++i){
			kp[(0+i)+(n+i)*ldkp] = ky[i]*Epsilon_inv[0]*kx[i];
			kp[(n+i)+(n+i)*ldkp] = omega2 - kx[i]*Epsilon_inv[0]*kx[i];
		}
	}else{
		RNP::TBLAS::CopyMatrix<'A'>(n,n, Epsilon_inv, n, &kp[0+0*ldkp], ldkp);
		RNP::TBLAS::CopyMatrix<'A'>(n,n, Epsilon_inv, n, &kp[n+0*ldkp], ldkp);
		RNP::TBLAS::CopyMatrix<'A'>(n,n, Epsilon_inv, n, &kp[0+n*ldkp], ldkp);
		RNP::TBLAS::CopyMatrix<'A'>(n,n, Epsilon_inv, n, &kp[n+n*ldkp], ldkp);

		for(size_t i = 0; i < n; ++i){
			RNP::TBLAS::Scale(n2, -ky[i], &kp[0+i*n2], 1);
		}
		for(size_t i = 0; i < n; ++i){
			RNP::TBLAS::Scale(n2, kx[i], &kp[0+(i+n)*n2], 1);
		}
		for(size_t i = 0; i < n; ++i){
			RNP::TBLAS::Scale(n2, ky[i], &kp[i+0*n2], ldkp);
		}
		for(size_t i = 0; i < n; ++i){
			RNP::TBLAS::Scale(n2, -kx[i], &kp[i+n+0*n2], ldkp);
		}
		const std::complex<double> omega2 = omega*omega;
		for(size_t i = 0; i < n2; ++i){
			kp[i+i*ldkp] += omega2;
		}
	}
}
static void MakeKPMatrix_real(
	double omega,
	size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv,
	double *kp_existing,
	double *kp,
	const size_t ldkp // leading dimension of kp, >= 2*n
){
	const size_t n2 = 2*n;
	if(NULL != kp_existing){
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, kp_existing,n2, kp,ldkp);
		return;
	}

	const double omega2 = omega*omega;

	for(size_t j = 0; j < n; ++j){
		for(size_t i = 0; i < n; ++i){
			kp[i+j*ldkp] = Epsilon_inv[i+j*n].real();
		}
	}

	RNP::TBLAS::CopyMatrix<'A'>(n,n, kp, ldkp, &kp[n+0*ldkp], ldkp);
	RNP::TBLAS::CopyMatrix<'A'>(n,n, kp, ldkp, &kp[0+n*ldkp], ldkp);
	RNP::TBLAS::CopyMatrix<'A'>(n,n, kp, ldkp, &kp[n+n*ldkp], ldkp);

	for(size_t i = 0; i < n; ++i){
		RNP::TBLAS::Scale(n2, -ky[i], &kp[0+i*n2], 1);
	}
	for(size_t i = 0; i < n; ++i){
		RNP::TBLAS::Scale(n2, kx[i], &kp[0+(i+n)*n2], 1);
	}
	for(size_t i = 0; i < n; ++i){
		RNP::TBLAS::Scale(n2, ky[i], &kp[i+0*n2], ldkp);
	}
	for(size_t i = 0; i < n; ++i){
		RNP::TBLAS::Scale(n2, -kx[i], &kp[i+n+0*n2], ldkp);
	}
	for(size_t i = 0; i < n2; ++i){
		kp[i+i*ldkp] += omega2;
	}
}

// kp = k-parallel matrix
// kp = omega^2 - Kappa = omega^2 - [  ky*epsinv*ky -ky*epsinv*kx ]
//                                  [ -kx*epsinv*ky  kx*epsinv*kx ]
//    = omega^2 + [ -ky ] [ epsinv ] [ ky -kx ]
//                [  kx ]

void MultKPMatrix(
	const char *trans,
	const std::complex<double> omega,
	const size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	const int Epsilon2_type,
	const std::complex<double> *kp, // kp if exists. Can be NULL
	const size_t ncols,
	const std::complex<double> *a,
	const size_t lda,
	std::complex<double> *y, // same size as a, y <- kp*a
	const size_t ldy
){
	const size_t n2 = 2*n;

	if(NULL != kp){
		if('N' == trans[0]){
			if(ncols > 1){
				RNP::TBLAS::MultMM<'N','N'>(n2,ncols,n2, std::complex<double>(1.),kp,n2, a,lda, std::complex<double>(0.),y,ldy);
			}else{
				RNP::TBLAS::MultMV<'N'>(n2,n2, std::complex<double>(1.),kp,n2, a,1, std::complex<double>(0.),y,1);
			}
		}else{
			if(ncols > 1){
				RNP::TBLAS::MultMM<'C','N'>(n2,ncols,n2, std::complex<double>(1.),kp,n2, a,lda, std::complex<double>(0.),y,ldy);
			}else{
				RNP::TBLAS::MultMV<'C'>(n2,n2, std::complex<double>(1.),kp,n2, a,1, std::complex<double>(0.),y,1);
			}
		}
		return;
	}

	const std::complex<double> omega2 = omega*omega;
	if(EPSILON2_TYPE_BLKDIAG1_SCALAR == Epsilon2_type || EPSILON2_TYPE_BLKDIAG2_SCALAR == Epsilon2_type){
		const std::complex<double> epsinv(
			'N' == trans[0] ? Epsilon_inv[0] : std::conj(Epsilon_inv[0])
		);
		for(size_t j = 0; j < ncols; ++j){
			for(size_t i = 0; i < n; ++i){
				std::complex<double> c = a[(0+i)+j*lda];
				std::complex<double> d = a[(n+i)+j*lda];
				y[(0+i)+j*ldy] = (omega2 - ky[i]*epsinv*ky[i]) * c + (         ky[i]*epsinv*kx[i]) * d;
				y[(n+i)+j*ldy] = (         kx[i]*epsinv*ky[i]) * c + (omega2 - kx[i]*epsinv*kx[i]) * d;
			}
		}
	}else{
		for(size_t j = 0; j < ncols; ++j){
			for(size_t i = 0; i < n; ++i){
				y[i+j*ldy] = ky[i]*a[(0+i)+j*lda] - kx[i]*a[(n+i)+j*lda];
			}
		}
		if('N' == trans[0]){
			if(ncols > 1){
				RNP::TBLAS::MultMM<'N','N'>(n,ncols,n, 1.,Epsilon_inv,n, &y[0+0*ldy],ldy, 0.,&y[n+0*ldy],ldy);
			}else{
				RNP::TBLAS::MultMV<'N'>(n,n, 1.,Epsilon_inv,n, &y[0+0*ldy],1, 0.,&y[n+0*ldy],1);
			}
		}else{
			if(ncols > 1){
				RNP::TBLAS::MultMM<'C','N'>(n,ncols,n, 1.,Epsilon_inv,n, &y[0+0*ldy],ldy, 0.,&y[n+0*ldy],ldy);
			}else{
				RNP::TBLAS::MultMV<'C'>(n,n, 1.,Epsilon_inv,n, &y[0+0*ldy],1, 0.,&y[n+0*ldy],1);
			}
		}

		for(size_t j = 0; j < ncols; ++j){
			for(size_t i = 0; i < n; ++i){
				y[(0+i)+j*ldy] = omega2 * a[(0+i)+j*lda] - ky[i]*y[n+i+j*ldy];
				y[(n+i)+j*ldy] = omega2 * a[(n+i)+j*lda] + kx[i]*y[n+i+j*ldy];
			}
		}
	}
}

void SolveLayerEigensystem_uniform(
	std::complex<double> omega,
	size_t n,
	const double *kx,
	const double *ky,
	std::complex<double> eps,
	std::complex<double> *q, // length 2*glist.n
	std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	std::complex<double> *phi // size (2*glist.n)^2
){
	const size_t n2 = 2*n;

	const std::complex<double> eta = 1./ eps;
	const std::complex<double> omega2 = omega*omega;

	// Fill kp
	if(NULL != kp){
		MakeKPMatrix(omega, n, kx, ky, &eta, EPSILON2_TYPE_BLKDIAG1_SCALAR, NULL, kp, n2);

#ifdef DUMP_MATRICES
	DUMP_STREAM << "kp:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,kp,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,kp,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
	}

	for(size_t i = 0; i < n; ++i){
		q[i] = eps*omega2 - kx[i]*kx[i] - ky[i]*ky[i];

		// Set the \hat{q} vector (diagonal matrix) while we're at it
		if(0 == omega.imag()){ // Not bandsolving
			if(std::abs(q[i].imag()) <= 4*DBL_EPSILON*std::abs(q[i].real())){
				if(q[i].real() >= 0){
					q[i] = std::sqrt(q[i].real());
				}else{
					q[i] = std::complex<double>(0, std::sqrt(-q[i].real()));
				}
			}else{
				q[i] = std::sqrt(q[i]);
				if(q[i].imag() < 0){
					q[i] = -q[i];
				}
			}
		}else{ // performing some kind of bandsolving, need to choose the appropriate branch
			if(q[i].real() < 0){
				// branch cut should be just below positive real axis
				q[i] = std::complex<double>(0,1) * std::sqrt(-q[i]);
			}else{
				// branch cut should be just below negative real axis
				// This is the default behavior for sqrt(std::complex)
				q[i] = std::sqrt(q[i]);
			}
		}

		/*
		if(std::complex<double>(0) == q[i]){
			q[i] = DBL_EPSILON*omega;
		}
		*/

		q[i+n] = q[i];
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "q:" << std::endl;
	RNP::IO::PrintVector(n2,q,1, DUMP_STREAM) << std::endl << std::endl;
#endif
	if(NULL != phi){
		RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., phi,n2);
#ifdef DUMP_MATRICES
	DUMP_STREAM << "phi:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,phi,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,phi,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
	}
}

#ifdef HAVE_LAPACK
void SolveLayerEigensystem_real(
	double omega,
	size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	const std::complex<double> *Epsilon2, // size (2*glist.n)^2 (dielectric/normal-field matrix)
	std::complex<double> *q, // length 2*glist.n
	std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	std::complex<double> *phi, // size (2*glist.n)^2
	std::complex<double> *work_,
	double *rwork_,
	size_t lwork
){
	const size_t n2 = 2*n;

	if((size_t)-1 == lwork){
		double tmp;
		RNP::Eigensystem_real(n2, NULL, n2, q, NULL, 1, phi, n2, &tmp, lwork);
		work_[0] = tmp + n2*n2;
		return;
	}else if(0 == lwork){
		lwork = n2*n2+2*n2;
	}

	double *work = (double*)work_;
	size_t eigenlwork;
	if(NULL == work_ || lwork < n2*n2+4*n2){
		lwork = (size_t)-1;
		double tmp;
		RNP::Eigensystem_real(n2, NULL, n2, q, NULL, 1, phi, n2, &tmp, lwork);
		eigenlwork = (size_t)tmp;
		work = (double*)rcwa_malloc(sizeof(double)*(eigenlwork + n2*n2));
	}else{
		eigenlwork = lwork - n2*n2;
	}
	double *op = work;
	double *eigenwork = op + n2*n2;

#ifdef DUMP_MATRICES
	DUMP_STREAM << "kx:" << std::endl;
	RNP::IO::PrintVector(n,kx,1, DUMP_STREAM) << std::endl;
	DUMP_STREAM << "ky:" << std::endl;
	RNP::IO::PrintVector(n,ky,1, DUMP_STREAM) << std::endl << std::endl;
#endif

#ifdef DUMP_MATRICES
	DUMP_STREAM << "Epsilon_inv:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n,n,Epsilon_inv,n, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n,Epsilon_inv,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	// Fill kp
	double *kp_use = (double*)(NULL != kp ? kp : phi);

	MakeKPMatrix_real(omega, n, kx, ky, Epsilon_inv, NULL, kp_use, n2);
#ifdef DUMP_MATRICES
	DUMP_STREAM << "kp:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,kp_use,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,&kp_use[0+(n2-1)*n2],1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif


#ifdef DUMP_MATRICES
	DUMP_STREAM << "Epsilon2:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,Epsilon2,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,Epsilon2,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	// Make the eigenoperator Epsilon2*kp - [kxkx, kxky; kykx, kyky]
	for(size_t j = 0; j < n2; ++j){
		for(size_t i = 0; i < n2; ++i){
	//		double sum = 0;
	//		for(size_t k = 0; k < n2; ++k){
	//			sum += Epsilon2[i+k*n2]*kp_use[k+j*n2];
	//		}
			op[i+j*n2] = RNP::TBLAS::Dot(n2, (double*)(&Epsilon2[i]), 2*n2, &kp_use[j*n2], 1);
		}
	}
	//RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., op,n2);
	//RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,Epsilon2,n2, kp_use,n2, 0.,op,n2);

	for(size_t i = 0; i < n; ++i){
		op[i+i*n2] -= kx[i]*kx[i];
	}
	for(size_t i = 0; i < n; ++i){
		op[i+n+i*n2] -= ky[i]*kx[i];
	}
	for(size_t i = 0; i < n; ++i){
		op[i+(i+n)*n2] -= kx[i]*ky[i];
	}
	for(size_t i = 0; i < n; ++i){
		op[i+n+(i+n)*n2] -= ky[i]*ky[i];
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "op:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,op,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintMatrix(n2,n2,op,n2, DUMP_STREAM) << std::endl << std::endl;
	//RNP::IO::PrintVector(n2,&op[0+(n2-1)*n2],1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

#ifdef DUMP_MATRICES
# ifdef DUMP_MATRICES_LARGE
	std::complex<double> *op_save = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>) * 2*n2*n2);
	std::complex<double> *op_temp = op_save + n2*n2;
	for(size_t j = 0; j < n2; ++j){
		for(size_t i = 0; i < n2; ++i){
			op_save[i+j*n2] = op[i+j*n2];
		}
	}
# endif
#endif
	int info = RNP::Eigensystem_real(n2, op, n2, q, NULL, 1, phi, n2, eigenwork, eigenlwork);
	if(0 != info){
		fprintf(stderr, "Layer eigensystem returned info = %d\n", info);
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "eigen info = " << info << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::TBLAS::Fill(n2*n2, 0., op_temp,1);
	RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,op_save,n2, phi,n2, 0.,op_temp,n2);
	for(size_t i = 0; i < n2; ++i){
		RNP::TBLAS::Axpy(n2, -q[i], &phi[0+i*n2],1, &op_temp[0+i*n2],1);
	}
	DUMP_STREAM << "eigensystem residual:" << std::endl;
	RNP::IO::PrintMatrix(n2,n2,op_temp,n2, DUMP_STREAM) << std::endl << std::endl;
	{
		for(size_t j = 0; j < n2; ++j){
			double sum = 0;
			for(size_t i = 0; i < n2; ++i){
				sum += fabs(op_temp[i+j*n2].real()) + fabs(op_temp[i+j*n2].imag());
			}
			DUMP_STREAM << " residual sum " << j << ": " << sum << std::endl;
		}
	}
	rcwa_free(op_save);
# endif
#endif

	for(size_t i = 0; i < n2; ++i){
		// Set the \hat{q} vector (diagonal matrix) while we're at it
		q[i] = std::sqrt(q[i]);
		if(q[i].imag() < 0){
			q[i] = -q[i];
		}
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "q:" << std::endl;
	RNP::IO::PrintVector(n2,q,1, DUMP_STREAM) << std::endl << std::endl;
	DUMP_STREAM << "phi:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,phi,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,phi,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	if(NULL == work_ || lwork < n2*n2+4*n2){
		rcwa_free(work);
	}
	
	// Re-make kp in complex arithmetic for future use
	std::complex<double> *kp_use_actual = (NULL != kp ? kp : phi);
	MakeKPMatrix(omega, n, kx, ky, Epsilon_inv, EPSILON2_TYPE_FULL, NULL, kp_use_actual, n2);
}
#endif // HAVE_LAPACK

void SolveLayerEigensystem(
	std::complex<double> omega,
	size_t n,
	const double *kx,
	const double *ky,
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	const std::complex<double> *Epsilon2, // size (2*glist.n)^2 (dielectric/normal-field matrix)
	int epstype,
	std::complex<double> *q, // length 2*glist.n
	std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	std::complex<double> *phi, // size (2*glist.n)^2
	std::complex<double> *work_,
	double *rwork_,
	size_t lwork
){
	const size_t n2 = 2*n;
	
#ifdef HAVE_LAPACK
	bool isreal = (n > 10);
	for(size_t j = 0; j < n && isreal; ++j){
		for(size_t i = 0; i < n; ++i){
			if(0 != Epsilon_inv[i+j*n].imag()){
				isreal = false;
				break;
			}
		}
	}
	for(size_t j = 0; j < n2 && isreal; ++j){
		for(size_t i = 0; i < n2; ++i){
			if(0 != Epsilon2[i+j*n2].imag()){
				isreal = false;
				break;
			}
		}
	}
	if(isreal && 0 == omega.imag() && EPSILON2_TYPE_FULL == epstype){
		SolveLayerEigensystem_real(
			omega.real(), n, kx, ky, Epsilon_inv, Epsilon2,
			q, kp, phi, work_, rwork_, lwork
		);
		return;
	}
#endif

	if((size_t)-1 == lwork){
		double dum;
		RNP::Eigensystem(n2, NULL, n2, q, NULL, 1, phi, n2, work_, &dum, lwork);
		work_[0] += n2*n2;
		return;
	}else if(0 == lwork){
		lwork = n2*n2+2*n2;
	}

	std::complex<double> *work = work_;
	size_t eigenlwork;
	if(NULL == work_ || lwork < n2*n2+2*n2){
		lwork = (size_t)-1;
		RNP::Eigensystem(n2, NULL, n2, q, NULL, 1, phi, n2, q, NULL, lwork);
		eigenlwork = (size_t)q[0].real();
		work = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*(eigenlwork + n2*n2));
	}else{
		eigenlwork = lwork - n2*n2;
	}
	std::complex<double> *op = work;
	std::complex<double> *eigenwork = op + n2*n2;

	double *rwork = rwork_;
	if(NULL == rwork_){
		rwork = (double*)rcwa_malloc(sizeof(double)*2*n2);
	}

#ifdef DUMP_MATRICES
	DUMP_STREAM << "kx:" << std::endl;
	RNP::IO::PrintVector(n,kx,1, DUMP_STREAM) << std::endl;
	DUMP_STREAM << "ky:" << std::endl;
	RNP::IO::PrintVector(n,ky,1, DUMP_STREAM) << std::endl << std::endl;
#endif

#ifdef DUMP_MATRICES
	DUMP_STREAM << "Epsilon_inv:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n,n,Epsilon_inv,n, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n,Epsilon_inv,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	// Fill kp
	std::complex<double> *kp_use = (NULL != kp ? kp : phi);

	MakeKPMatrix(omega, n, kx, ky, Epsilon_inv, epstype, NULL, kp_use, n2);
#ifdef DUMP_MATRICES
	DUMP_STREAM << "kp:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,kp_use,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,&kp_use[0+(n2-1)*n2],1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif


#ifdef DUMP_MATRICES
	DUMP_STREAM << "Epsilon2:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,Epsilon2,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,Epsilon2,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	// Make the eigenoperator Epsilon2*kp - [kxkx, kxky; kykx, kyky]
	RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., op,n2);
	RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,Epsilon2,n2, kp_use,n2, 0.,op,n2);

	for(size_t i = 0; i < n; ++i){
		op[i+i*n2] -= kx[i]*kx[i];
	}
	for(size_t i = 0; i < n; ++i){
		op[i+n+i*n2] -= ky[i]*kx[i];
	}
	for(size_t i = 0; i < n; ++i){
		op[i+(i+n)*n2] -= kx[i]*ky[i];
	}
	for(size_t i = 0; i < n; ++i){
		op[i+n+(i+n)*n2] -= ky[i]*ky[i];
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "op:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,op,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintMatrix(n2,n2,op,n2, DUMP_STREAM) << std::endl << std::endl;
	//RNP::IO::PrintVector(n2,&op[0+(n2-1)*n2],1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

#ifdef DUMP_MATRICES
# ifdef DUMP_MATRICES_LARGE
	std::complex<double> *op_save = (std::complex<double> *)rcwa_malloc(sizeof(std::complex<double>) * 2*n2*n2);
	std::complex<double> *op_temp = op_save + n2*n2;
	RNP::TBLAS::CopyMatrix<'A'>(n2,n2, op,n2, op_save,n2);
# endif
#endif
	int info = RNP::Eigensystem(n2, op, n2, q, NULL, 1, phi, n2, eigenwork, rwork, eigenlwork);
	if(0 != info){
		fprintf(stderr, "Layer eigensystem returned info = %d\n", info);
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "eigen info = " << info << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::TBLAS::Fill(n2*n2, 0., op_temp,1);
	RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,op_save,n2, phi,n2, 0.,op_temp,n2);
	for(size_t i = 0; i < n2; ++i){
		RNP::TBLAS::Axpy(n2, -q[i], &phi[0+i*n2],1, &op_temp[0+i*n2],1);
	}
	DUMP_STREAM << "eigensystem residual:" << std::endl;
	RNP::IO::PrintMatrix(n2,n2,op_temp,n2, DUMP_STREAM) << std::endl << std::endl;
	rcwa_free(op_save);
# endif
#endif

	for(size_t i = 0; i < n2; ++i){
		// Set the \hat{q} vector (diagonal matrix) while we're at it
		if(0 == omega.imag()){ // Not bandsolving
			q[i] = std::sqrt(q[i]);
			if(q[i].imag() < 0){
				q[i] = -q[i];
			}
		}else{ // performing some kind of bandsolving, need to choose the appropriate branch
			if(q[i].real() < 0){
				// branch cut should be just below positive real axis
				q[i] = std::complex<double>(0,1) * std::sqrt(-q[i]);
			}else{
				// branch cut should be just below negative real axis
				// This is the default behavior for sqrt(std::complex)
				q[i] = std::sqrt(q[i]);
			}
		}
	}
#ifdef DUMP_MATRICES
	DUMP_STREAM << "q:" << std::endl;
	RNP::IO::PrintVector(n2,q,1, DUMP_STREAM) << std::endl << std::endl;
	DUMP_STREAM << "phi:" << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n2,n2,phi,n2, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n2,phi,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	if(NULL == work_ || lwork < n2*n2+2*n2){
		rcwa_free(work);
	}
	if(NULL == rwork_){
		rcwa_free(rwork);
	}
}


void InitSMatrix(
	size_t n,
	std::complex<double> *S // size (4*n)^2
){
	const size_t n4 = 4*n;
	RNP::TBLAS::SetMatrix<'A'>(n4,n4, 0.,1., S, n4);
}
void GetSMatrix(
	size_t nlayers,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const double *thickness, // list of thicknesses
	const std::complex<double> **q, // list of q vectors
	const std::complex<double> **Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int *epstype,
	const std::complex<double> **kp,
	const std::complex<double> **phi,
	std::complex<double> *S, // size (4*n)^2
	std::complex<double> *work_,
	size_t *iwork,
	size_t lwork
){
	if(0 == nlayers){ return; }
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;

	if((size_t)-1 == lwork){
		work_[0] = n4*(n4+1);
		return;
	}
	std::complex<double> *work = work_;
	if(NULL == work_ || lwork < n4*(n4+1)){
		work = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*(n4*(n4+1)));
	}
	size_t *pivots = iwork;
	if(NULL == iwork){
		pivots = (size_t*)rcwa_malloc(sizeof(size_t)*n4);
	}

	RNP::TBLAS::SetMatrix<'A'>(n4,n4, 0.,1., S, n4);

	std::complex<double> *t1 = work;
	std::complex<double> *t2 = t1 + n2*n2;
	std::complex<double> *in1 = t2 + n2*n2;
	std::complex<double> *in2 = in1 + n2*n2;
	std::complex<double> *d1 = in2 + n2*n2;
	std::complex<double> *d2 = d1 + n2;

	for(size_t l = 0; l < nlayers-1; ++l){
		size_t lp1 = l+1;
		if(lp1 >= nlayers){ lp1 = l; }

		// Make the interface matrices
		if((lp1 == l) || (q[l] == q[lp1] && ((NULL != kp[l] && kp[l] == kp[lp1]) || Epsilon_inv[l] == Epsilon_inv[lp1]) && phi[l] == phi[lp1])){
			// This is a trivial interface, set to identity
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., in1, n2);
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., in2, n2);
		}else{
			// The interface matrix is the inverse of the mode-to-field matrix of layer l
			// times the mode-to-field matrix of layer l+1 (lp1).
			// The mode-to-field matrix is of the form
			// [ B -B ] where A = phi
			// [ A  A ] where B = kp*phi*inv(diag(q)) = G*A/q
			// So we want
			// 0.5 * [  iBl  iAl ] [ Blp1 -Blp1 ]
			//       [ -iBl  iAl ] [ Alp1  Alp1 ]
			// Multiplying out gives
			// 0.5 * [ P+Q P-Q ] // where P = iAl*Alp1, and i in front means inverse
			//       [ P-Q P+Q ] // where Q = iBl*Blp1
			// Making P is easy, since A is a single matrix.
			// Making Q is as follows:
			// Q = iBl*Blp1
			//   = ql*iAl*iGl * Gl*Alp1*iqlp1
			// We will only store I11 and I21
			/*
			{
				std::complex<double> *Ml = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*n4*n4);
				std::complex<double> *Mlp1 = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*n4*n4);

				RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., &Ml[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, &Ml[n2+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, &Ml[n2+n2*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, &Ml[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, 1./q[l][i], &Ml[0+(i+n2)*n4], 1);
				}
				RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,kp[l],n2, &Ml[0+n2*n4],n4, 0.,&Ml[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Ml[0+0*n4],n4, &Ml[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, -1., &Ml[0+(i+n2)*n4], 1);
				}

				RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., &Mlp1[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, &Mlp1[n2+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, &Mlp1[n2+n2*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, &Mlp1[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, 1./q[lp1][i], &Mlp1[0+(i+n2)*n4], 1);
				}
				RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,kp[lp1],n2, &Mlp1[0+n2*n4],n4, 0.,&Mlp1[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Mlp1[0+0*n4],n4, &Mlp1[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, -1., &Mlp1[0+(i+n2)*n4], 1);
				}

				//RNP::LinearSolve<'N'>(n4, n4, Ml, n4, Mlp1, n4, NULL, pivots);
#ifdef DUMP_MATRICES
			DUMP_STREAM << "M(" << l << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
			RNP::IO::PrintMatrix(n4,n4,Ml,n4, DUMP_STREAM) << std::endl << std::endl;
# else
			RNP::IO::PrintVector(n4,Ml,1, DUMP_STREAM) << std::endl << std::endl;
# endif
			DUMP_STREAM << "M(" << lp1 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
			RNP::IO::PrintMatrix(n4,n4,Mlp1,n4, DUMP_STREAM) << std::endl << std::endl;
# else
			RNP::IO::PrintVector(n4,Mlp1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

				rcwa_free(Mlp1);
				rcwa_free(Ml);
			}
			*/
			// Make Bl in t1
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., t1,n2);
			{
				if(NULL == phi[l]){
					if(NULL == kp[l]){
						MakeKPMatrix(omega, n, kx, ky, Epsilon_inv[l], epstype[l], kp[l], t1,n2);
					}else{
						RNP::TBLAS::CopyMatrix<'A'>(n2,n2, kp[l],n2, t1,n2);
					}
				}else{
					MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv[l], epstype[l], kp[l], n2, phi[l],n2, t1,n2);
				}
				//for(size_t i = 0; i < n2; ++i){
				//	RNP::TBLAS::Scale(n2, 1./q[l][i], &t1[0+i*n2], 1);
				//}
			}
#ifdef DUMP_MATRICES
			DUMP_STREAM << "Bl(" << l << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
			RNP::IO::PrintMatrix(n2,n2,t1,n2, DUMP_STREAM) << std::endl << std::endl;
# else
			RNP::IO::PrintVector(n2,t1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
			// Make Blp1 in in1
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., in1,n2);
			{
				if(NULL == phi[lp1]){
					if(NULL == kp[lp1]){
						MakeKPMatrix(omega, n, kx, ky, Epsilon_inv[lp1], epstype[lp1], kp[lp1], in1,n2);
					}else{
						RNP::TBLAS::CopyMatrix<'A'>(n2,n2, kp[lp1],n2, in1,n2);
					}
				}else{
					MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv[lp1], epstype[lp1], kp[lp1], n2, phi[lp1],n2, in1,n2);
				}
				//for(size_t i = 0; i < n2; ++i){
				//	RNP::TBLAS::Scale(n2, 1./q[lp1][i], &in1[0+i*n2], 1);
				//}
			}
#ifdef DUMP_MATRICES
		DUMP_STREAM << "Bl(" << l+1 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n2,n2,in1,n2, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n2,in1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
			int solve_info;
			// Make Q in in1
			//RNP::LinearSolve<'N'>(n2, n2, t1, n2, in1, n2, &solve_info, pivots);
			SingularLinearSolve(n2,n2,n2, t1,n2, in1,n2, DBL_EPSILON);
			// Now perform the diagonal scalings
			for(size_t i = 0; i < n2; ++i){
				RNP::TBLAS::Scale(n2, q[l][i], &in1[i+0*n2], n2);
			}
			{
				double maxel = 0;
				for(size_t i = 0; i < n2; ++i){
					double el = std::abs(q[lp1][i]);
					if(el > maxel){ maxel = el; }
				}
				for(size_t i = 0; i < n2; ++i){
					double el = std::abs(q[lp1][i]);
					if(el < DBL_EPSILON * maxel){
						RNP::TBLAS::Scale(n2, 0., &in1[0+i*n2], 1);
					}else{
						RNP::TBLAS::Scale(n2, 1./q[lp1][i], &in1[0+i*n2], 1);
					}
				}
			}

			// Make P in in2
			if(NULL == phi[lp1]){
				RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., in2,n2);
			}else{
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, in2,n2);
			}
			if(NULL != phi[l]){
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, t1,n2);
				RNP::LinearSolve<'N'>(n2, n2, t1, n2, in2, n2, &solve_info, pivots);
			}

			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, in2,n2, t1,n2); // in2 = P, t1 = P, in1 = Q
			RNP::TBLAS::Axpy(n2*n2, -1., in1,1, in2,1); // in2 = P-Q, t1 = P, in1 = Q
			RNP::TBLAS::Axpy(n2*n2, 1., t1,1, in1,1); // in2 = P+Q, t1 = P, in1 = P+Q
			RNP::TBLAS::Scale(n2*n2, 0.5, in1,1);
			RNP::TBLAS::Scale(n2*n2, 0.5, in2,1);
		}
#ifdef DUMP_MATRICES
		DUMP_STREAM << "Interface1(" << l+1 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n2,n2,in1,n2, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n2,in1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
		DUMP_STREAM << "Interface2(" << l+1 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n2,n2,in2,n2, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n2,in2,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

		for(size_t i = 0; i < n2; ++i){
			d1[i] = std::exp(q[l  ][i] * std::complex<double>(0,thickness[l  ]));
			d2[i] = std::exp(q[lp1][i] * std::complex<double>(0,thickness[lp1]));
		}

		// Make S11
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, -1.,&S[0+n2*n4],n4, in2,n2, 0.,t1,n2); // t1 = -S12 I21
		for(size_t i = 0; i < n2; ++i){ // t1 = -f_l S12 I21
			RNP::TBLAS::Scale(n2, d1[i], &t1[i+0*n2], n2);
		}
		RNP::TBLAS::Axpy(n2*n2, 1., in1,1, t1,1); // t1 = (I11 - f_l S12 I21)

		RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., t2,n2);
		int solve_info;
		RNP::LinearSolve<'N'>(n2, n2, t1, n2, t2, n2, &solve_info, pivots); // t2 = (I11 - f_l S12 I21)^{-1}

		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &S[0+0*n4],n4, t1,n2);
		for(size_t i = 0; i < n2; ++i){ // t1 = f_l S11
			RNP::TBLAS::Scale(n2, d1[i], &t1[i+0*n2], n2);
		}
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,t2,n2, t1,n2, 0.,&S[0+0*n4],n4);
		// S11 is done, and we need to hold on to t2 = (I11 - f_l S12 I21)^{-1}

		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,&S[0+n2*n4],n4, in1,n2, 0.,t1,n2); // t1 = S12 I22
		for(size_t i = 0; i < n2; ++i){ // t1 = f_l S12 I22
			RNP::TBLAS::Scale(n2, d1[i], &t1[i+0*n2], n2);
		}
		RNP::TBLAS::Axpy(n2*n2, -1., in2,1, t1,1); // t1 = f_l S12 I22 - I12
		for(size_t i = 0; i < n2; ++i){ // t1 = (f_l S12 I22 - I12) f_{l+1}
			RNP::TBLAS::Scale(n2, d2[i], &t1[0+i*n2], 1);
		}
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,t2,n2, t1,n2, 0.,&S[0+n2*n4],n4);
		// S12 done, and t2 can be reused

		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,&S[n2+n2*n4],n4, in2,n2, 0.,t1,n2); // t1 = S22 I21
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,t1,n2, &S[0+0*n4],n4, 1.,&S[n2+0*n4],n4);
		// S21 done, need to keep t1 = S22 I21

		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,&S[n2+n2*n4],n4, in1,n2, 0.,t2,n2); // t2 = S22 I22
		for(size_t i = 0; i < n2; ++i){ // t2 = S22 I22 f_{l+1}
			RNP::TBLAS::Scale(n2, d2[i], &t2[0+i*n2], 1);
		}
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, t2,n2, &S[n2+n2*n4],n4);
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,t1,n2, &S[0+n2*n4],n4, 1.,&S[n2+n2*n4],n4);

#ifdef DUMP_MATRICES
		DUMP_STREAM << "S(1," << l+2 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n4,n4,S,n4, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n4,S,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
	}
	if(NULL == work_ || lwork < n4*(n4+1)){
		rcwa_free(work);
	}
	if(NULL == iwork){
		rcwa_free(pivots);
	}
}


int SolveAll(
	size_t nlayers,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const double *thickness, // list of thicknesses
	const std::complex<double> **q, // list of q vectors
	const std::complex<double> **Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int *epstype,
	const std::complex<double> **kp,
	const std::complex<double> **phi,
	std::complex<double> *ab, // length 2*n*nlayers
	std::complex<double> *work_, // length lwork
	size_t *iwork, // length n2*nlayers
	size_t lwork // set to -1 for query into iwork[0], at least 6*n2^2*nlayers
){
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;
	const size_t n22 = n2*n2;
	const size_t minwork = 6*n22*nlayers;
	if((size_t)-1 == lwork){
		iwork[0] = minwork;
		return 0;
	}
	typedef std::complex<double> doublecomplex;
	doublecomplex *work = work_;
	if(0 == lwork && NULL == work_){
		work = (doublecomplex*)rcwa_malloc(sizeof(doublecomplex) * minwork);
		if(NULL == work){
printf("Could not allocate %d\n", (int)minwork); fflush(stdout);
			return 1;
		}
	}else if(lwork < minwork){
		// error: not enough workspace
		return -15;
	}
	
	// Set up temporaries
	doublecomplex *t1 = work;
	doublecomplex *t2 = t1 + n22;
	doublecomplex *in1 = t2 + n22;
	doublecomplex *in2 = in1 + n22;
	
	for(size_t ip = 0, iq = 1; iq < nlayers; ip=iq++){
		// Set pointers to current row of matrices
		doublecomplex *row = work+6*n22*iq;
		
		// Precondition for current loop iteration:
		//   Qprev from previous iteration was formed correctly
		
		/*
		// Assemble interface T-matrix
		//   Assemble mode-to-field operator for p in row
		if(0 == p){
			// TODO: Generate mode-to-field operator for p
			
			// Factor and invert P
			if(layer p is uniform){
			}else if(layer p is z-uncoupled){
			}else{
				LUFactor(n4, n4, row, n4, ipiv);
				LUInvert(n4, n4, row, n4, ipiv);
			}
		}else{ // Mode-to-field operator for p was generated and factored in previous iteration
			// Only invert P (in-place)
			if(layer p is uniform){
			}else if(layer p is z-uncoupled){
			}else{
				LUInvert(n4, n4, row, n4, ipiv);
			}
		}
		//   Assemble mode-to-field operator for q in rownext
		// TODO: Generate mode-to-field operator for q
		if(layer q is uniform){
		}else if(layer q is z-uncoupled){
		}else{
			LUFactor(n4, n4, rownext, n4, ipiv);
		}
		// In-place multiply q into p
		LUMultRight(n4, n4, rownext, n4, ipiv);
		*/
		
		// Make the interface matrices
		if((ip == iq) || (q[ip] == q[iq] && ((NULL != kp[ip] && kp[ip] == kp[iq]) || Epsilon_inv[ip] == Epsilon_inv[iq]) && phi[ip] == phi[iq])){
			// This is a trivial interface, set to identity
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., in1, n2);
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., in2, n2);
		}else{
			// The interface matrix is the inverse of the mode-to-field matrix of layer l
			// times the mode-to-field matrix of layer l+1 (lp1).
			// The mode-to-field matrix is of the form
			// [ B -B ] where A = phi
			// [ A  A ] where B = kp*phi*inv(diag(q)) = G*A/q
			// So we want
			// 0.5 * [  iBl  iAl ] [ Blp1 -Blp1 ]
			//       [ -iBl  iAl ] [ Alp1  Alp1 ]
			// Multiplying out gives
			// 0.5 * [ P+Q P-Q ] // where P = iAl*Alp1, and i in front means inverse
			//       [ P-Q P+Q ] // where Q = iBl*Blp1
			// Making P is easy, since A is a single matrix.
			// Making Q is as follows:
			// Q = iBl*Blp1
			//   = ql*iAl*iGl * Gl*Alp1*iqlp1
			// We will only store I11 and I21
			/*
			{
				std::complex<double> *Ml = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*n4*n4);
				std::complex<double> *Mlp1 = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*n4*n4);

				RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., &Ml[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, &Ml[n2+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, &Ml[n2+n2*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[l],n2, &Ml[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, 1./q[l][i], &Ml[0+(i+n2)*n4], 1);
				}
				RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,kp[l],n2, &Ml[0+n2*n4],n4, 0.,&Ml[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Ml[0+0*n4],n4, &Ml[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, -1., &Ml[0+(i+n2)*n4], 1);
				}

				RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., &Mlp1[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, &Mlp1[n2+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, &Mlp1[n2+n2*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[lp1],n2, &Mlp1[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, 1./q[lp1][i], &Mlp1[0+(i+n2)*n4], 1);
				}
				RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, 1.,kp[lp1],n2, &Mlp1[0+n2*n4],n4, 0.,&Mlp1[0+0*n4],n4);
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Mlp1[0+0*n4],n4, &Mlp1[0+n2*n4],n4);
				for(size_t i = 0; i < n2; ++i){
					RNP::TBLAS::Scale(n2, -1., &Mlp1[0+(i+n2)*n4], 1);
				}

				//RNP::LinearSolve<'N'>(n4, n4, Ml, n4, Mlp1, n4, NULL, pivots);
#ifdef DUMP_MATRICES
			DUMP_STREAM << "M(" << l << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
			RNP::IO::PrintMatrix(n4,n4,Ml,n4, DUMP_STREAM) << std::endl << std::endl;
# else
			RNP::IO::PrintVector(n4,Ml,1, DUMP_STREAM) << std::endl << std::endl;
# endif
			DUMP_STREAM << "M(" << lp1 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
			RNP::IO::PrintMatrix(n4,n4,Mlp1,n4, DUMP_STREAM) << std::endl << std::endl;
# else
			RNP::IO::PrintVector(n4,Mlp1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

				rcwa_free(Mlp1);
				rcwa_free(Ml);
			}
			*/
			
			//// Begin copy from GetSMatrix
			// Make Bl in t1
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., t1,n2);
			{
				if(NULL == phi[ip]){
					if(NULL == kp[ip]){
						MakeKPMatrix(omega, n, kx, ky, Epsilon_inv[ip], epstype[ip], kp[ip], t1,n2);
					}else{
						RNP::TBLAS::CopyMatrix<'A'>(n2,n2, kp[ip],n2, t1,n2);
					}
				}else{
					MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv[ip], epstype[ip], kp[ip], n2, phi[ip],n2, t1,n2);
				}
			}
#ifdef DUMP_MATRICES
			DUMP_STREAM << "Bl(" << ip << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
			RNP::IO::PrintMatrix(n2,n2,t1,n2, DUMP_STREAM) << std::endl << std::endl;
# else
			RNP::IO::PrintVector(n2,t1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
			// Make Blp1 in in1
			RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,0., in1,n2);
			{
				if(NULL == phi[iq]){
					if(NULL == kp[iq]){
						MakeKPMatrix(omega, n, kx, ky, Epsilon_inv[iq], epstype[iq], kp[iq], in1,n2);
					}else{
						RNP::TBLAS::CopyMatrix<'A'>(n2,n2, kp[iq],n2, in1,n2);
					}
				}else{
					MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv[iq], epstype[iq], kp[iq], n2, phi[iq],n2, in1,n2);
				}
			}
#ifdef DUMP_MATRICES
		DUMP_STREAM << "Bl(" << iq << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n2,n2,in1,n2, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n2,in1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif
			int solve_info;
			// Make Q in in1
			RNP::LinearSolve<'N'>(n2, n2, t1, n2, in1, n2, &solve_info, iwork);
			//SingularLinearSolve(n2,n2,n2, t1,n2, in1,n2, DBL_EPSILON);
			// Now perform the diagonal scalings
			for(size_t i = 0; i < n2; ++i){
				RNP::TBLAS::Scale(n2, q[ip][i], &in1[i+0*n2], n2);
			}
			{
				double maxel = 0;
				for(size_t i = 0; i < n2; ++i){
					double el = std::abs(q[iq][i]);
					if(el > maxel){ maxel = el; }
				}
				for(size_t i = 0; i < n2; ++i){
					double el = std::abs(q[iq][i]);
					if(el < DBL_EPSILON * maxel){
						RNP::TBLAS::Scale(n2, 0., &in1[0+i*n2], 1);
					}else{
						RNP::TBLAS::Scale(n2, 1./q[iq][i], &in1[0+i*n2], 1);
					}
				}
			}

			// Make P in in2
			if(NULL == phi[iq]){
				RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., in2,n2);
			}else{
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[iq],n2, in2,n2);
			}
			if(NULL != phi[ip]){
				RNP::TBLAS::CopyMatrix<'A'>(n2,n2, phi[ip],n2, t1,n2);
				RNP::LinearSolve<'N'>(n2, n2, t1, n2, in2, n2, &solve_info, iwork);
			}

			RNP::TBLAS::CopyMatrix<'A'>(n2,n2, in2,n2, t1,n2); // in2 = P, t1 = P, in1 = Q
			RNP::TBLAS::Axpy(n2*n2, -1., in1,1, in2,1); // in2 = P-Q, t1 = P, in1 = Q
			RNP::TBLAS::Axpy(n2*n2, 1., t1,1, in1,1); // in2 = P+Q, t1 = P, in1 = P+Q
			RNP::TBLAS::Scale(n2*n2, 0.5, in1,1);
			RNP::TBLAS::Scale(n2*n2, 0.5, in2,1);
		}
#ifdef DUMP_MATRICES
		DUMP_STREAM << "Interface1(" << iq << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n2,n2,in1,n2, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n2,in1,1, DUMP_STREAM) << std::endl << std::endl;
# endif
		DUMP_STREAM << "Interface2(" << iq << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
		RNP::IO::PrintMatrix(n2,n2,in2,n2, DUMP_STREAM) << std::endl << std::endl;
# else
		RNP::IO::PrintVector(n2,in2,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

		//PrintMatrix("in1", n2, n2, in1, n2);
		//PrintMatrix("in2", n2, n2, in2, n2);

		doublecomplex *Sba = row;
		doublecomplex *Saa = Sba+n2;
		doublecomplex *Sbb = row+n2*n4;
		doublecomplex *Sab = Sbb+n2;
		
//printf("Preparing to exchange T to S\n"); fflush(stdout);

		// Exchange to interface S-matrix, and also swap block rows
		Copy(n2,n2, in1,n2, Saa,n4);
		Invert(n2, Saa,n4, Sbb, n22, iwork); // Use Sbb as workspace
		Mult(n2,  1., in2,n2, Saa,n4, 0., Sba,n4);
		Mult(n2, -1., Saa,n4, in2,n2, 0., Sab,n4);
		Copy(n2,n2, in1,n2, Sbb,n4);
		Mult(n2,  1., in2,n2, Sab,n4, 1., Sbb,n4);

//printf("Exchanged T to S\n"); fflush(stdout);

//		PrintMatrix("Saa", n2, n2, Saa, n4);
//		PrintMatrix("Sab", n2, n2, Sab, n4);
//		PrintMatrix("Sba", n2, n2, Sba, n4);
//		PrintMatrix("Sbb", n2, n2, Sbb, n4);
		
		// Fold in propagation phases into the S-matrix
		for(size_t j = 0; j < n2; ++j){
			doublecomplex sp = std::exp(q[ip][j] * std::complex<double>(0,thickness[ip]));
			doublecomplex sq = std::exp(q[iq][j] * std::complex<double>(0,thickness[iq]));
//printf("Scale factor ip j = %f + %f i\n", sp.real(), sp.imag());
//printf("Scale factor iq j = %f + %f i\n", sq.real(), sq.imag());
			Scale(n4, sp, &row[0+j*n4]);
			Scale(n4, sq, &row[0+(j+n2)*n4]);
		}

//		PrintMatrix("Saa with phase", n2, n2, Saa, n4);
//		PrintMatrix("Sab with phase", n2, n2, Sab, n4);
//		PrintMatrix("Sba with phase", n2, n2, Sba, n4);
//		PrintMatrix("Sbb with phase", n2, n2, Sbb, n4);

//printf("Applied phases\n"); fflush(stdout);
	}
	///////////// End copy from GetSMatrix
	
	// At this point, the workspace is divided into portions of length 6*n22
	// for each layer. The interface scattering matrix between layers i and i+1
	// begins at offset 6*n22*(i+1) in the workspace. The contents of each block
	// has the form:
	//   [ Sba Sbb ]
	//   [ Saa Sab ]
	// The system of equations has instead a block row of the form:
	//   [ -Sba  I   0  -Sbb ]
	//   [ -Saa  0   I  -Sab ]
	// Assembling all such rows, one obtains a rectangular matrix of size
	// 4*n*(m-1) by 4*n*m where m is the number of layers.
	// The labeling of the rows is [ b0, a1, b1, a2, b2, ... am ]
	// The labeling of the columns is [ a0, b0, a1, b1, a2, b2, ... am, bm ]
	// The first and last block columns correspond to boundary conditions and
	// hence must be moved to the RHS. This produces in the end a square matrix.
	
	// Prepare the RHS
	Mult(n4, n2, 1., &work[6*n22], n4, &ab[0], 0, t1);
	Copy(n4, t1, 1, &ab[n2], 1);
	Mult(n4, n2, 1., &work[6*n22*(nlayers-1)+n2*n4], n4, &ab[(nlayers-1)*n4+n2], 0, t1);
	Copy(n4, t1, 1, &ab[(nlayers-1)*n4-n2], 1);
	
//	PrintMatrix("RHS0", n4,nlayers, ab, n4);
	
	// At this point, we can compute the quantities involved in the LDU decomposition
	// of the system matrix. The initial diagonal pivot block is
	//   [ I    0  | -Sbb  0 ]
	//   [ 0    I  | -Sab  0 ]
	//   [ --------+-------- ]
	//   [ 0  -Sba |         ]
	//   [ 0  -Saa |         ]
	// Notice that the upper and lower diagonal blocks come from different interfaces.
	// The first LDU factorization step is trivial, and there are m-1 in total.
	// Therefore, we need only iterate through the m-2 inner layers.
	
	for(size_t i = 1; i < nlayers; i++){
//printf("i = %d\n", (int)i);
		// Set pointers to current row of matrices
		doublecomplex *row = work+6*n22*i;
		doublecomplex *P = row + 4*n22;
		doublecomplex *Q = P + n22;
		doublecomplex *Pprev = P - 6*n22;
		doublecomplex *Qprev = Pprev + n22;
			
		const doublecomplex *Sbb = row-6*n22+n2*n4;
		const doublecomplex *Sab = Sbb+n2;
		
		const doublecomplex *Sba = row;
		const doublecomplex *Saa = Sba+n2;
		size_t *ipivP = iwork + n2*i;
		
		// Perform LDU factorization step
		if(1 == i){ // First step is trivial
			SetMatrix(n2, n2, 1., 0., P, n2);
			ZeroMatrix(n2, n2, Q, n2);
		}else{
			Copy(n2, n2, Sab, n4, t1, n2); // t1 = Sab
			Copy(n2, n2, Sbb, n4, Q, n2); // Q = Sbb, t1 = Sab
			LUSolve(n2, n2, Pprev, n2, ipivP-n2, Q, n2); // Q = inv(P) Sbb, t1 = Sab
//	PrintMatrix("t1(Sab)", n2,n2, Sab, n4);
//	PrintMatrix("         Sbb", n2,n2, Sbb, n4);
//	PrintMatrix("invPprev Sbb", n2,n2, Q, n2);
//	PrintMatrix("Qprev", n2,n2, Qprev, n2);
			Mult(n2, 1., Qprev, n2, Q, n2, -1., t1, n2); // t1 = Q*inv(P)*Sbb - Sab
//	PrintMatrix("Qprev invPprev Sbb", n2,n2, t1, n2);
			Mult(n2, 1., Saa, n4, t1, n2, 0., Q, n2);
	//PrintMatrix("Saa", n2,n2, Saa, n4);
	//PrintMatrix("Saa Qprev invPprev Sbb", n2,n2, Q, n2);
			
			Mult(n2, 1., Sba, n4, t1, n2, 0., P, n2);
			for(size_t j = 0; j < n2; ++j){
				P[j+j*n2] += 1.;
			}
			
	//		PrintMatrix("Q", n2,n2, Q, n2);
	//		PrintMatrix("P", n2,n2, P, n2);
		}
//	PrintMatrix("P", n2,n2, P, n2);
		LUFactor(n2, P, n2, ipivP);
//	PrintMatrix("Pfactored", n2,n2, P, n2);
//	for(size_t i = 0; i < n2; ++i){
//		printf(" %d", ((int*)ipivP)[i]);
//	} printf("\n");

//printf("LDU step completed\n"); fflush(stdout);
	}
	
	// Now use LDU factorization to solve
	//   L D U x = b
	//   x = U \ (D \ (L \ b))
	//
	// Forward pass for D \ (L \ b)
	for(size_t j = 1; j+1 < nlayers; ++j){
		// Set pointers to current row of matrices
		doublecomplex *row = work+6*n22*j;
		doublecomplex *P = row + 4*n22;
		doublecomplex *Q = P + n22;
		doublecomplex *Pprev = P - 6*n22;
		doublecomplex *Qprev = Pprev + n22;
			
		doublecomplex *Sbb = row+n2*n4;
		doublecomplex *Sab = Sbb+n2;
		
		doublecomplex *Sba = row+6*n22;
		doublecomplex *Saa = Sba+n2;
		
		doublecomplex *aj = &ab[j*n4];
		doublecomplex *bjm1 = aj-n2;
		doublecomplex *bjp1 = aj+n2;
		size_t *ipivP = iwork + n2*j;
		
		// sub-diagonal block:
		//   [     P               |  ]
		//   [     Q            I  |  ]
		//   [ --------------------+- ]
		//   [ Sba Q inv(P)   -Sba |  ]
		//   [ Saa Q inv(P)   -Saa |  ]
		Copy(n4, bjm1, 1, t1, 1);
//PrintMatrix("bjm1", n4,1, bjm1, n4);
//PrintMatrix("Applied matrix", n4,n2, Sba, n4);
		Mult(n4, n2, 1., Sba, n4, &t1[n2], 1., bjp1);
		if(j > 1){ // Second iteration should use the fact that Q = 0
			LUSolve(n2, 1, P, n2, ipivP, t1, n2);
			Mult(n2, n2, -1., Q, n2, t1, 0., &t1[n2]);
			Mult(n4, n2,  1., Sba, n4, &t1[n2], 1., bjp1);
		}
	}
//printf("Forward pass completed\n"); fflush(stdout);
//PrintMatrix("RHS1", n4,nlayers, ab, n4);
	
	// Diagonal pass
	for(size_t j = 2; j < nlayers; ++j){
		// Set pointers to current row of matrices
		const doublecomplex *row = work+6*n22*j;
		const doublecomplex *P = row + 4*n22;
		const doublecomplex *Q = P + n22;
		size_t *ipivP = iwork + n2*j;
		doublecomplex *aj = &ab[j*n4];
		doublecomplex *bjm1 = aj-n2;
		// Diaonal block:
		//   [ P     0 ]
		//   [ Q     I ]
		// Inverse:
		//   [    inv(P)   0 ] [ bjm1 ]
		//   [ -Q inv(P)   I ] [ aj   ]
		LUSolve(n2, 1, P, n2, ipivP, bjm1, n2);
//PrintMatrix("bjm1", n2,1, bjm1, n2);
		Mult(n2, n2, -1., Q, n2, bjm1, 1., aj);
//PrintMatrix("Q", n2,n2, Q, n2);
	}
//printf("Diagonal pass completed\n"); fflush(stdout);
//PrintMatrix("RHS2", n4,nlayers, ab, n4);
	
	// Backward pass for U \ ...
	for(size_t j = nlayers-2; j > 0; --j){
		// Set pointers to current row of matrices
		doublecomplex *row = work+6*n22*j;
		doublecomplex *P = row + 4*n22;
		doublecomplex *Q = P + n22;
		doublecomplex *Pprev = P - 6*n22;
		doublecomplex *Qprev = Pprev + n22;
			
		doublecomplex *Sbb = row+n2*n4;
		doublecomplex *Sab = Sbb+n2;
		
		doublecomplex *Sba = row+6*n22;
		doublecomplex *Saa = Sba+n2;
		
		doublecomplex *aj = &ab[j*n4];
		doublecomplex *bjm1 = aj-n2;
		doublecomplex *bjp1 = aj+n2;
		size_t *ipivP = iwork + n2*j;
		
		// super-diagonal block:
		//   [ P   0 | -inv(P) Sbb            0 ]
		//   [ Q   I | Q inv(P) Sbb - Sab     0 ]
		//   [ ------+------------------------- ]
		//   [       |                          ]
		Mult(n2, n2, 1., Sbb, n4, bjp1, 0., t1);
		if(j > 0){
			LUSolve(n2, 1, P, n2, ipivP, t1, n2);
		}
		Axpy(n2, 1., t1, 1, bjm1, 1);
		if(j > 0){
			Mult(n2, n2, -1., Q, n2, t1, 1., aj);
		}
		Mult(n2, n2,  1., Sab, n4, bjp1, 1., aj);
	}
	
//printf("Backward pass completed\n"); fflush(stdout);
//PrintMatrix("RHS3", n4,nlayers, ab, n4);
	
	if(NULL == work_){
		rcwa_free(work);
	}
	return 0;
}

int SolveInterior(
	size_t nlayers,
	size_t which_layer,
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const double *thickness, // list of thicknesses
	const std::complex<double> **q, // list of q vectors
	const std::complex<double> **Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int *epstype,
	const std::complex<double> **kp,
	const std::complex<double> **phi,
	std::complex<double> *a0, // length 2*n
	std::complex<double> *bN, // length 2*n
	std::complex<double> *ab, // length 4*n
	std::complex<double> *work_, // length lwork
	size_t *iwork, // length n2
	size_t lwork // set to -1 for query into work[0], at least 2*(4*n)^2 + 2*(2*n) + 4*n*(4*n+1)
){
	if(0 == nlayers){ return-1; }
	if(which_layer >= nlayers){ return -2; }

	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;
	const size_t lwork_GetSMatrix = 4*n*(4*n+1);
	const size_t lwork_needed = 2*n4*n4 + 2*n2 + lwork_GetSMatrix;

	if((size_t)-1 == lwork){
		work_[0] = lwork_needed;
		return 0;
	}
	std::complex<double> *work = work_;
	if(NULL == work_ || lwork < lwork_needed){
		work = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*lwork_needed);
	}
	size_t *pivots = iwork;
	if(NULL == iwork){
		pivots = (size_t*)rcwa_malloc(sizeof(size_t)*n2);
	}

	std::complex<double> *S0l = work;
	std::complex<double> *SlN = S0l + n4*n4;
	std::complex<double> *S11a0 = SlN + n4*n4;
	std::complex<double> *S22bN = S11a0 + n2;
	std::complex<double> *work_GetSMatrix = S22bN + n2;
	std::complex<double> *al = ab;
	std::complex<double> *bl = al+n2;
	std::complex<double> *temp = S0l;
	size_t ldtemp = n4;

	int info;

	GetSMatrix(which_layer+1, n, kx, ky, omega,
		thickness, q, Epsilon_inv, epstype, kp, phi,
		S0l, work_GetSMatrix, pivots, lwork_GetSMatrix);
	GetSMatrix(nlayers-which_layer, n, kx, ky, omega,
		thickness+which_layer, q+which_layer, Epsilon_inv+which_layer, epstype+which_layer, kp+which_layer, phi+which_layer,
		SlN, work_GetSMatrix, pivots, lwork_GetSMatrix);

#ifdef DUMP_MATRICES
	if(NULL != a0){
		DUMP_STREAM << "a0 = " << std::endl;
		RNP::IO::PrintVector(n2,a0,1, DUMP_STREAM) << std::endl << std::endl;
	}
	if(NULL != bN){
		DUMP_STREAM << "bN = " << std::endl;
		RNP::IO::PrintVector(n2,bN,1, DUMP_STREAM) << std::endl << std::endl;
	}
#endif
#ifdef DUMP_MATRICES
	DUMP_STREAM << "S0l(0," << which_layer << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n4,n4,S0l,n4, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n4,S0l,1, DUMP_STREAM) << std::endl << std::endl;
# endif
	DUMP_STREAM << "SlN(" << which_layer << "," << nlayers-1 << ") = " << std::endl;
# ifdef DUMP_MATRICES_LARGE
	RNP::IO::PrintMatrix(n4,n4,SlN,n4, DUMP_STREAM) << std::endl << std::endl;
# else
	RNP::IO::PrintVector(n4,SlN,1, DUMP_STREAM) << std::endl << std::endl;
# endif
#endif

	// both solutions only depend on the products S11(0,l)*a0 and S22(l,N)*bN
	if(NULL != a0){
		RNP::TBLAS::MultMV<'N'>(n2, n2, std::complex<double>(1.0), &S0l[0+0*n4], n4,
			a0, 1,
			std::complex<double>(0.0), S11a0, 1);
	}else{
		RNP::TBLAS::Fill(n2, 0., S11a0, 1);
	}
	if(NULL != bN){
		RNP::TBLAS::MultMV<'N'>(n2, n2, std::complex<double>(1.0), &SlN[n2+n2*n4], n4,
			bN, 1,
			std::complex<double>(0.0), S22bN, 1);
	}else{
		RNP::TBLAS::Fill(n2, 0., S22bN, 1);
	}

	// We overwrite the upper left submatrix S11(0,l) since it's not needed anymore
	// temp is set to S0l for this reason.

	// Compute -S_12(0,l)S_21(l,N)
	RNP::TBLAS::MultMM<'N','N'>(n2, n2, n2, std::complex<double>(-1.0), &S0l[0+n2*n4], n4,
		&SlN[n2+0*n4], n4,
		std::complex<double>(0.0), temp, ldtemp);
	for(size_t i = 0; i < n2; ++i){
		temp[i+i*ldtemp] += 1.;
	} // temp = (1 - S_12(0,l)S_21(l,N))

	RNP::TBLAS::MultMV<'N'>(n2, n2, std::complex<double>(1.0), &S0l[0+n2*n4], n4,
		S22bN, 1,
		std::complex<double>(0.0), al, 1); // al = S_12(0,l)S_22(l,N)bN
	RNP::TBLAS::Axpy(n2, std::complex<double>(1.0), S11a0, 1, al, 1); // al = S_11(0,l)*a0 + S_12(0,l)S_22(l,N)bN

	RNP::LinearSolve<'N'>(n2, 1, temp, ldtemp, al, n2, &info, pivots);
	// al done

	// Make the other matrix
	// Compute S_21(l,N)S_12(0,l)
	RNP::TBLAS::MultMM<'N','N'>(n2, n2, n2, std::complex<double>(-1.0), &SlN[n2+0*n4], n4,
		&S0l[0+n2*n4], n4,
		std::complex<double>(0.0), temp, ldtemp);
	for(size_t i = 0; i < n2; ++i){
		temp[i+i*ldtemp] += 1.;
	} // temp = (1 - S_21(l,N)S_12(0,l))

	RNP::TBLAS::MultMV<'N'>(n2, n2, std::complex<double>(1.0), &SlN[n2+0*n4], n4,
		S11a0, 1,
		std::complex<double>(0.0), bl, 1); // bl = S_21(l,N)S_11(0,l)a0
	RNP::TBLAS::Axpy(n2, std::complex<double>(1.0), S22bN, 1, bl, 1); // bl = S_21(l,N)S_11(0,l)a0 + S_22(l,N)bN
	RNP::LinearSolve<'N'>(n2, 1, temp, ldtemp, bl, n2, &info, pivots);

#ifdef DUMP_MATRICES
	DUMP_STREAM << "al:" << std::endl;
	RNP::IO::PrintVector(n2,al,1, DUMP_STREAM) << std::endl << std::endl;
	DUMP_STREAM << "bl:" << std::endl;
	RNP::IO::PrintVector(n2,bl,1, DUMP_STREAM) << std::endl << std::endl;
#endif

	if(NULL == work_ || lwork < lwork_needed){
		rcwa_free(work);
	}
	if(NULL == iwork){
		rcwa_free(pivots);
	}
	return 0;
}



void TranslateAmplitudes(
	size_t n, // glist.n
	const std::complex<double> *q, // length 2*glist.n
	double thickness,
	double dz,
	std::complex<double> *ab
){
	const size_t n2 = 2*n;
	for(size_t i = 0; i < n2; ++i){
		std::complex<double> iq = std::complex<double>(0,1)*q[i];
		// When extrapolating, the exponentials can be huge, so we need to check
		// against overflow. This doesn't solve the problem completely, but it helps.
		if(0. != ab[i]){
			ab[i]    *= std::exp(iq*dz);
		}
		if(0. != ab[i+n2]){
			ab[i+n2] *= std::exp(iq*(thickness-dz));
		}
	}
}

void GetZPoyntingFlux(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int epstype,
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *forward,
	std::complex<double> *backward,
	std::complex<double> *work
){
	const size_t n2 = 2*n;

	std::complex<double> *a2 = work;
	if(NULL == work){
		a2 = (std::complex<double> *)rcwa_malloc(sizeof(std::complex<double>) * 4*n2);
	}
	std::complex<double> *b2 = a2 + n2;
	std::complex<double> *a3 = b2 + n2;
	std::complex<double> *b3 = a3 + n2;

	memcpy(a2, ab, sizeof(std::complex<double>) * 2*n2);
	for(size_t i = 0; i < n2; ++i){ a2[i] /= (omega*q[i]); }
	for(size_t i = 0; i < n2; ++i){ b2[i] /= (omega*q[i]); }
	if(NULL == phi){
		RNP::TBLAS::CopyMatrix<'A'>(n2,2, a2,n2, a3,n2);
	}else{
		RNP::TBLAS::MultMM<'N','N'>(n2,2,n2, std::complex<double>(1.),phi,n2, a2,n2, std::complex<double>(0.), a3,n2);
	}
	MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv, epstype, kp, 2, a3,n2, a2,n2);
	if(NULL == phi){
		RNP::TBLAS::CopyMatrix<'A'>(n2,2, a2,n2, a3,n2);
	}else{
		RNP::TBLAS::MultMM<'C','N'>(n2,2,n2, std::complex<double>(1.),phi,n2, a2,n2, std::complex<double>(0.), a3,n2);
	}

	// At this point a3[n2] is alpha and b3[n2] is beta,
	// ab[n2] is a and (ab+n2)[n2] is b

	*forward  =  RNP::TBLAS::ConjugateDot(n2, ab     ,1, a3,1).real();
	*backward = -RNP::TBLAS::ConjugateDot(n2, &ab[n2],1, b3,1).real();
	std::complex<double> diff = 0.5*(RNP::TBLAS::ConjugateDot(n2, &ab[n2],1, a3,1) - RNP::TBLAS::ConjugateDot(n2, b3,1, ab,1));
	*forward  += diff;
	*backward += std::conj(diff);

	if(NULL == work){
		rcwa_free(a2);
	}
}

void GetZPoyntingFluxComponents(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int epstype,
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *forward,
	std::complex<double> *backward,
	std::complex<double> *work
){
	const size_t n2 = 2*n;

	std::complex<double> *a2 = work;
	if(NULL == work){
		a2 = (std::complex<double> *)rcwa_malloc(sizeof(std::complex<double>) * 4*n2);
	}
	std::complex<double> *b2 = a2 + n2;
	std::complex<double> *a3 = b2 + n2;
	std::complex<double> *b3 = a3 + n2;

	memcpy(a2, ab, sizeof(std::complex<double>) * 2*n2);
	for(size_t i = 0; i < n2; ++i){ a2[i] /= (omega*q[i]); }
	for(size_t i = 0; i < n2; ++i){ b2[i] /= (omega*q[i]); }
	if(NULL == phi){
		RNP::TBLAS::CopyMatrix<'A'>(n2,2, a2,n2, a3,n2);
	}else{
		RNP::TBLAS::MultMM<'N','N'>(n2,2,n2, std::complex<double>(1.),phi,n2, a2,n2, std::complex<double>(0.), a3,n2);
	}
	MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv, epstype, kp, 2, a3,n2, a2,n2);
	// At this point, a2[n2] contains alpha_e, b2[n2] contains beta_e

	if(NULL == phi){
		RNP::TBLAS::CopyMatrix<'A'>(n2,2, ab,n2, a3,n2);
	}else{
		RNP::TBLAS::MultMM<'N','N'>(n2,2,n2, std::complex<double>(1.),phi,n2, ab,n2, std::complex<double>(0.), a3,n2);
	}
	// At this point, a3[n2] contains alpha_h, b3[n2] contains beta_h

	for(size_t i = 0; i < n; ++i){
		forward[i] = 0;
		backward[i] = 0;
		for(size_t j = 0; j < n2; j+=n){
			const size_t k = i+j;
			forward[i]  += std::real(std::conj(a2[k])*a3[k]);
			backward[i] -= std::real(std::conj(b2[k])*b3[k]);
			std::complex<double> diff = 0.5*(std::conj(b3[k])*a2[k] - std::conj(b2[k])*a3[k]);
			forward[i]  += diff;
			backward[i] += std::conj(diff);
		}
	}

	if(NULL == work){
		rcwa_free(a2);
	}
}

static void GetInPlaneFieldVector(
	size_t n, // glist.n
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *Epsilon_inv, // size (glist.n)^2; inv of usual dielectric Fourier coupling matrix
	int epstype,
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *eh // length 8*2*glist.n
){
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;

	for(size_t i = 0; i < n2; ++i){
		eh[4*n2+i] = ab[i]/(omega*q[i]);
		eh[5*n2+i] = -ab[i+n2]/(omega*q[i]);
	}
	RNP::TBLAS::Copy(n4, ab,1, &eh[6*n2],1);
	if(NULL == phi){
		RNP::TBLAS::CopyMatrix<'A'>(n2,4, &eh[4*n2],n2, &eh[0],n2);
	}else{
		RNP::TBLAS::MultMM<'N','N'>(n2,4,n2, std::complex<double>(1.),phi,n2, &eh[4*n2],n2, std::complex<double>(0.),&eh[0],n2);
	}
	// At this point, the 4 columns of &eh[0] (length n2) are
	// [ phi.inv(q).a   -phi.inv(q).b   phi.a   phi.b ]
	MultKPMatrix("N", omega, n, kx, ky, Epsilon_inv, epstype, kp, 2, &eh[0],n2, &eh[4*n2],n2);
	RNP::TBLAS::Axpy(n2, std::complex<double>(1.),&eh[2*n2],1, &eh[3*n2],1);
	RNP::TBLAS::Axpy(n2, std::complex<double>(1.),&eh[5*n2],1, &eh[4*n2],1);
	/*
	const std::complex<double> *hx  = &eh[3*n2+0];
	const std::complex<double> *hy  = &eh[3*n2+n];
	const std::complex<double> *ney = &eh[4*n2+0];
	const std::complex<double> *ex  = &eh[4*n2+n];
	*/
}

void GetFieldAtPoint(
	size_t n, // glist.n
	const double *kx,
	const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2, non NULL for efield != NULL
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	const double r[2], // coordinates within layer
	std::complex<double> efield[3],
	std::complex<double> hfield[3],
	std::complex<double> *work // 8*n2
){
	const std::complex<double> z_zero(0.);
	const std::complex<double> z_one(1.);
	const size_t n2 = 2*n;

	std::complex<double> *eh = work;
	if(NULL == work){
		eh = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>) * 8*n2);
	}

	GetInPlaneFieldVector(n, kx, ky, omega, q, epsilon_inv, epstype, kp, phi, ab, eh);
	const std::complex<double> *hx  = &eh[3*n2+0];
	const std::complex<double> *hy  = &eh[3*n2+n];
	const std::complex<double> *ney = &eh[4*n2+0];
	const std::complex<double> *ex  = &eh[4*n2+n];

	std::complex<double> fE[3], fH[3];
	fE[0] = 0;
	fE[1] = 0;
	fE[2] = 0;
	fH[0] = 0;
	fH[1] = 0;
	fH[2] = 0;

	if(NULL != efield && NULL != epsilon_inv){
		for(size_t i = 0; i < n; ++i){
			eh[i] = (ky[i]*hx[i] - kx[i]*hy[i]);
		}
		if(EPSILON2_TYPE_BLKDIAG1_SCALAR == epstype || EPSILON2_TYPE_BLKDIAG2_SCALAR == epstype){
			RNP::TBLAS::Scale(n, epsilon_inv[0], eh,1);
			RNP::TBLAS::Copy(n, eh,1, &eh[n], 1);
		}else{
			RNP::TBLAS::MultMV<'N'>(n,n, z_one,epsilon_inv,n, eh,1, z_zero,&eh[n],1);
		}
	}

	for(size_t i = 0; i < n; ++i){
		const double theta = (kx[i]*r[0] + ky[i]*r[1]);

		const std::complex<double> phase(cos(theta),sin(theta));
		fH[0] += hx[i]*phase;
		fH[1] += hy[i]*phase;
		fE[0] += ex[i]*phase;
		fE[1] -= ney[i]*phase;
		fH[2] += (kx[i] * -ney[i] - ky[i] * ex[i]) * (phase / omega);
		fE[2] += eh[n+i] * (phase / omega);
	}

	if(NULL != efield && NULL != epsilon_inv){
		efield[0] = fE[0];
		efield[1] = fE[1];
		efield[2] = fE[2];
	}
	if(NULL != hfield){
		hfield[0] = fH[0];
		hfield[1] = fH[1];
		hfield[2] = fH[2];
	}

	if(NULL == work){
		rcwa_free(eh);
	}
}
void GetFieldOnGrid(
	size_t n, // glist.n
	int *G,
	const double *kx, const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2, non NULL for efield != NULL || kp == NULL
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	const size_t nxy[2], // number of points per lattice direction
	const double *xy0,
	std::complex<double> *efield,
	std::complex<double> *hfield
){
	const std::complex<double> z_zero(0.);
	const std::complex<double> z_one(1.);
	const size_t n2 = 2*n;
	const size_t N = nxy[0]*nxy[1];
	const int nxyoff[2] = { (int)(nxy[0]/2), (int)(nxy[1]/2) };
	int inxy[2] = { (int)nxy[0], (int)nxy[1] };
	int inxy_rev[2] = { (int)nxy[1], (int)nxy[0] };

	std::complex<double> *eh = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>) * 8*n2);

	GetInPlaneFieldVector(n, kx, ky, omega, q, epsilon_inv, epstype, kp, phi, ab, eh);
	const std::complex<double> *hx  = &eh[3*n2+0];
	const std::complex<double> *hy  = &eh[3*n2+n];
	const std::complex<double> *ney = &eh[4*n2+0];
	const std::complex<double> *ex  = &eh[4*n2+n];

	std::complex<double> *from[6];
	std::complex<double> *to[6];
	fft_plan plan[6];
	for(unsigned i = 0; i < 6; ++i){
		from[i] = fft_alloc_complex(N);
		to[i] = fft_alloc_complex(N);
		memset(from[i], 0, sizeof(std::complex<double>) * N);
		plan[i] = fft_plan_dft_2d(inxy_rev, from[i], to[i], 1);
	}

	for(size_t i = 0; i < n; ++i){
		eh[i] = (ky[i]*hx[i] - kx[i]*hy[i]);
	}
	if(EPSILON2_TYPE_BLKDIAG1_SCALAR == epstype || EPSILON2_TYPE_BLKDIAG2_SCALAR == epstype){
		RNP::TBLAS::Scale(n, epsilon_inv[0], eh,1);
		RNP::TBLAS::Copy(n, eh,1, &eh[n], 1);
	}else{
		RNP::TBLAS::MultMV<'N'>(n,n, z_one,epsilon_inv,n, eh,1, z_zero,&eh[n],1);
	}
	for(size_t i = 0; i < n; ++i){
		const int iu = G[2*i+0];
		const int iv = G[2*i+1];
		if(
			(nxyoff[0] - (int)nxy[0] < iu && iu <= nxyoff[0]) &&
			(nxyoff[1] - (int)nxy[1] < iv && iv <= nxyoff[1])
		){
			//todo: const std::complex<double> shift_phase(-i 2pi Lk.G.xy0);
			const int ii = (iu >= 0 ? iu : iu + nxy[0]);
			const int jj = (iv >= 0 ? iv : iv + nxy[1]);
			from[0][ii+jj*nxy[0]] = hx[i];
			from[1][ii+jj*nxy[0]] = hy[i];
			from[2][ii+jj*nxy[0]] = (kx[i] * ney[i] + ky[i] * ex[i]) / omega;
			from[3][ii+jj*nxy[0]] = ex[i];
			from[4][ii+jj*nxy[0]] = -ney[i];
			from[5][ii+jj*nxy[0]] = eh[n+i] / omega;
		}
	}

	for(unsigned i = 0; i < 6; ++i){
		fft_plan_exec(plan[i]);
	}

	for(size_t j = 0; j < nxy[1]; ++j){
		for(size_t i = 0; i < nxy[0]; ++i){
			hfield[3*(i+j*nxy[0])+0] = to[0][i+j*nxy[0]];
			hfield[3*(i+j*nxy[0])+1] = to[1][i+j*nxy[0]];
			hfield[3*(i+j*nxy[0])+2] = to[2][i+j*nxy[0]];
			efield[3*(i+j*nxy[0])+0] = to[3][i+j*nxy[0]];
			efield[3*(i+j*nxy[0])+1] = to[4][i+j*nxy[0]];
			efield[3*(i+j*nxy[0])+2] = to[5][i+j*nxy[0]];
		}
	}

	for(unsigned i = 0; i < 6; ++i){
		fft_plan_destroy(plan[i]);
		fft_free(to[i]);
		fft_free(from[i]);
	}
	rcwa_free(eh);
}

void GetZStressTensorIntegral(
	size_t n, // glist.n
	const double *kx,
	const double *ky,
	std::complex<double> omega,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2
	const std::complex<double> *epsilon2, // size (2*glist.n)^2
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> integral[3],
	std::complex<double> *work
){
	const size_t n2 = 2*n;
	const std::complex<double> z_zero(0.);
	const std::complex<double> z_one(1.);

	std::complex<double> *eh = work;
	if(NULL == work){
		eh = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>) * 8*n2);
	}

	GetInPlaneFieldVector(n, kx, ky, omega, q, epsilon_inv, epstype, kp, phi, ab, eh);
	std::complex<double> *ez  = &eh[1*n2+n];
	std::complex<double> *dz  = &eh[2*n2+0];
	std::complex<double> *hz  = &eh[2*n2+n];
	const std::complex<double> *hx  = &eh[3*n2+0];
	const std::complex<double> *hy  = &eh[3*n2+n];
	const std::complex<double> *ney = &eh[4*n2+0];
	const std::complex<double> *ex  = &eh[4*n2+n];
	std::complex<double> *ndy = &eh[5*n2+0];
	std::complex<double> *dx  = &eh[5*n2+n];

	integral[0] = 0;
	integral[1] = 0;
	integral[2] = 0;

	for(size_t i = 0; i < n; ++i){
		hz[i] = (kx[i] * -ney[i] - ky[i] * ex[i]) / omega;
		dz[i] = (ky[i]*hx[i] - kx[i]*hy[i]) / omega;
	}
	if(EPSILON2_TYPE_BLKDIAG1_SCALAR == epstype || EPSILON2_TYPE_BLKDIAG2_SCALAR == epstype){
		RNP::TBLAS::Copy(n, dz,1, ez, 1);
		RNP::TBLAS::Scale(n, epsilon_inv[0], ez,1);
	}else{
		RNP::TBLAS::MultMV<'N'>(n,n, z_one,epsilon_inv,n, dz,1, z_zero,ez,1);
	}
	RNP::TBLAS::MultMV<'N'>(n2,n2, z_one,epsilon2,n2, ney,1, z_zero,ndy,1);

	for(size_t i = 0; i < n; ++i){
		std::complex<double> cdz = std::conj(dz[i]);
		std::complex<double> cbz = std::conj(hz[i]);
		integral[0] += ex[i]*cdz;
		integral[1] +=-ney[i]*cdz;
		integral[2] += 0.5*ez[i]*cdz;
		integral[2] -= 0.5*ney[i]*std::conj(ndy[i]);
		integral[2] -= 0.5*ex[i]*std::conj(dx[i]);

		integral[0] += hx[i]*cbz;
		integral[1] += hy[i]*cbz;
		integral[2] += 0.5*hz[i]*std::conj(hz[i]);
		integral[2] -= 0.5*hy[i]*std::conj(hy[i]);
		integral[2] -= 0.5*hx[i]*std::conj(hx[i]);
	}

	if(NULL == work){
		rcwa_free(eh);
	}
}

std::complex<double> zsinc(const std::complex<double> &z){
	if(std::abs(z) < DBL_EPSILON){
		return 1.0;
	}else{
		return std::sin(z)/z;
	}
}

void GetLayerVolumeIntegral(
	char which,
	size_t n, // glist.n
	const double *kx,
	const double *ky,
	std::complex<double> omega,
	double thickness,
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2
	const std::complex<double> *epsilon2, // size (2*glist.n)^2
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	std::complex<double> *integral,
	std::complex<double> *work
){
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;
	const std::complex<double> omega2 = omega*omega;
	const double aomega2(std::norm(omega));
	const double omega_2(1./aomega2);
	std::complex<double> omega2_2(1.);
	if(0 != omega.imag()){
		omega2_2 = omega_2*omega2;
	}
	const std::complex<double> z_zero(0.);
	const std::complex<double> z_one(1.);

	std::complex<double> *Q = work;
	if(NULL == work){
		Q = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*n4*n4);
	}
	RNP::TBLAS::SetMatrix<'A'>(n4,n4, 0.,0., Q,n4);

	// Form CsC (in Q)
	if('E' == which || 'U' == which){
		RNP::TBLAS::CopyMatrix<'A'>(n2,n2, epsilon2,n2, &Q[0+0*n4],n4);
	}else if('e' == which){
		RNP::TBLAS::SetMatrix<'A'>(n2,n2, 0.,1., &Q[0+0*n4],n4);
	}
	if('H' == which || 'U' == which){
		for(size_t i = 0; i < n; ++i){
			Q[i+i*n4] += omega_2*kx[i]*kx[i];
		}
		for(size_t i = 0; i < n; ++i){
			Q[i+n+i*n4] += omega_2*ky[i]*kx[i];
		}
		for(size_t i = 0; i < n; ++i){
			Q[i+(i+n)*n4] += omega_2*kx[i]*ky[i];
		}
		for(size_t i = 0; i < n; ++i){
			Q[i+n+(i+n)*n4] += omega_2*ky[i]*ky[i];
		}
	}
	// At this point, top block of CsC is done
	if('E' == which || 'U' == which){
		MakeKPMatrix(omega, n, kx, ky, epsilon_inv, epstype, kp, &Q[n2+n2*n4],n4);
		for(size_t i = 0; i < n2; ++i){
			RNP::TBLAS::Scale(n2, -omega_2, &Q[n2+(i+n2)*n4],1);
			Q[i+n2+(i+n2)*n4] += z_one;
			Q[i+n2+(i+n2)*n4] *= omega2_2;
		}
	}else if('e' == which){
		if(EPSILON2_TYPE_BLKDIAG1_SCALAR == epstype || EPSILON2_TYPE_BLKDIAG2_SCALAR == epstype){
			RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,std::norm(epsilon_inv[0]), &Q[n2+n2*n4],n4);
		}else{
			RNP::TBLAS::MultMM<'C','N'>(n,n,n, z_one,epsilon_inv,n, epsilon_inv,n, z_zero,&Q[n2+n2*n4],n4);
		}
		RNP::TBLAS::CopyMatrix<'A'>(n,n, &Q[n2+n2*n4],n4, &Q[n2+n+n2*n4],n4);
		RNP::TBLAS::CopyMatrix<'A'>(n,n, &Q[n2+n2*n4],n4, &Q[n2+(n2+n)*n4],n4);
		RNP::TBLAS::CopyMatrix<'A'>(n,n, &Q[n2+n2*n4],n4, &Q[n2+n+(n2+n)*n4],n4);
		for(size_t i = 0; i < n; ++i){
			RNP::TBLAS::Scale(n2, ky[i], &Q[n2+0+(n2+i)*n4],1);
			RNP::TBLAS::Scale(n2, -kx[i], &Q[n2+0+(n2+n+i)*n4],1);
			RNP::TBLAS::Scale(n2, ky[i], &Q[n2+i+n2*n4],n4);
			RNP::TBLAS::Scale(n2, -kx[i], &Q[n2+n+i+n2*n4],n4);
		}
	}
	if('H' == which || 'U' == which){
		for(size_t i = 0; i < n2; ++i){
			Q[i+n2+(i+n2)*n4] += z_one;
		}
	}

	// Apply the amplitude-to-field matrix M to either side
	// CsC = [ A   ]   M = [ kp.phi.inv(q)   -kp.phi.inv(q) ] = [ fie -fie ]
	//       [   B ]       [     phi                phi     ]   [ phi  phi ]
	// Q = M^H * CsC * M = [  S+T -S+T ], S = fie^H * A * fie
	//                     [ -S+T  S+T ]  T = phi^H * B * phi
	// integral = int_0^thickness [x(z)^H * Q * x(z)] dz, x(z) = [a*exp(iqz);b*exp(iq(d-z))]

	// Make T
	if(NULL != phi){
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, z_one,&Q[n2+n2*n4],n4, phi,n2, z_zero,&Q[0+n2*n4],n4);
		RNP::TBLAS::MultMM<'C','N'>(n2,n2,n2, z_one,phi,n2, &Q[0+n2*n4],n4, z_zero, &Q[n2+n2*n4],n4);
	}
	// At this point: Q = [ A   ]
	//                    [   T ]
	if(NULL == kp){
		MakeKPMatrix(omega, n, kx, ky, epsilon_inv, epstype, kp, &Q[0+n2*n4],n4);
	}
	const std::complex<double> *kp_use = (kp != NULL ? kp : &Q[0+n2*n4]);
	const size_t ldkp = (kp != NULL ? n2 : n4);
	RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, z_one,&Q[0+0*n4],n4, kp_use,ldkp, z_zero,&Q[n2+0*n4],n4);
	RNP::TBLAS::MultMM<'C','N'>(n2,n2,n2, z_one,kp_use,ldkp, &Q[n2+0*n4],n4, z_zero, &Q[0+0*n4],n4);
	if(NULL != phi){
		RNP::TBLAS::MultMM<'N','N'>(n2,n2,n2, z_one,&Q[0+0*n4],n4, phi,n2, z_zero,&Q[n2+0*n4],n4);
		RNP::TBLAS::MultMM<'C','N'>(n2,n2,n2, z_one,phi,n2, &Q[n2+0*n4],n4, z_zero, &Q[0+0*n4],n4);
	}
	for(size_t i = 0; i < n2; ++i){
		RNP::TBLAS::Scale(n2, 1./(omega*q[i]), &Q[0+i*n4],1);
		RNP::TBLAS::Scale(n2, 1./(omega*q[i]), &Q[i+0*n4],n4);
	}
	// At this point: Q = [ S   ]
	//                    [   T ]
	RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Q[n2+n2*n4],n4, &Q[n2+0*n4],n4);
	RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Q[n2+n2*n4],n4, &Q[0+n2*n4],n4);
	for(size_t i = 0; i < n2; ++i){
		RNP::TBLAS::Axpy(n2, -1., &Q[0+i*n4],1, &Q[n2+i*n4],1);
		RNP::TBLAS::Axpy(n2, -1., &Q[0+i*n4],1, &Q[0+(i+n2)*n4],1);
	}
	// At this point: Q = [  S   T-S ]
	//                    [ T-S   T  ]
	for(size_t i = 0; i < n2; ++i){
		RNP::TBLAS::Axpy(n2, 1., &Q[0+i*n4],1, &Q[n2+(i+n2)*n4],1);
	}
	RNP::TBLAS::CopyMatrix<'A'>(n2,n2, &Q[n2+n2*n4],n4, &Q[0+0*n4],n4);
	// Q done

	// The integral we need to construct is
	// \int_0^d x(z)^H * Q * x(z) dz, where x(z) = [ a*e^{iqz}; b*e^{iq(d-z)} ]
	// We are summing terms like
	// \int_0^d x_i(z)^* x_j(z) Q_{ij} dz
	// The generic form is
	// Q_{ij} ab[i]^* ab[j] \int_0^d e^{-iq^* L_i(z)} e^{iq L_j(z)},
	// where L_i(z) is either z or d-z, and s_i is either 1 or -1.
	// Note that L_i(d) = d or 0, so L_i(d) = d*(1+s_i)/2, and L_i(0) = d*(1-s_i)/2
	// The integral simplifies to
	// \frac{1}{i(s_j q - s_i q^*)} [ e^{iq L_j(d)-iq^* L_i(d)} - e^{iq L_j(0)-iq^* L_i(0)} ]
	// \frac{1}{i(s_j q - s_i q^*)} [ e^{id/2 (q (s_j+1)-q^* (s_i+1))} - e^{id/2 (q (1-s_j)-q^* (1-s_i))} ]
	// \frac{1}{i(s_j q - s_i q^*)} [ e^{id/2 (q-q^*)} e^{id/2 (s_j q - s_i q^*)} - e^{id/2 (q-q^*)} e^{-id/2 (s_j q - s_i q^*)} ]
	// \frac{e^{id/2 (q-q^*)} }{i(s_j q - s_i q^*)} [ e^{id/2 (s_j q - s_i q^*)} - e^{-id/2 (s_j q - s_i q^*)} ]
	// \frac{e^{id/2 (q-q^*)} }{i(s_j q - s_i q^*)} 2 i sin[ d/2 (s_j q - s_i q^*) ]
	// \frac{d e^{id/2 (q-q^*)} }{d/2 (s_j q - s_i q^*)} sin[ d/2 (s_j q - s_i q^*) ]
	// d e^{id/2 (q-q^*)} sinc[ d/2 (s_j q - s_i q^*) ], where sinc(x) = sin(x)/x

	*integral = 0;
	const double d_2 = 0.5*thickness;
	for(size_t j = 0; j < n4; ++j){
		const double sj = (j < n2) ? 1 : -1;
		const size_t qj = (j < n2) ? j : j-n2;
		for(size_t i = 0; i < n4; ++i){
			const double si = (i < n2) ? 1 : -1;
			const size_t qi = (i < n2) ? i : i-n2;
			std::complex<double> term = Q[i+j*n4] * std::conj(ab[i]) * ab[j];
			term *= thickness*std::exp(std::complex<double>(0,d_2) * (q[qj] - std::conj(q[qi])));
			term *= zsinc(d_2 * (sj*q[qj] - si*std::conj(q[qi])));
			*integral += term;
		}
	}
	if(NULL == work){
		rcwa_free(Q);
	}
}


// returns Integral[ abs(XX)^2, z] for XX = ex, ey, ez, hx, hy, hz
void GetLayerZIntegral(
	size_t n, // glist.n
	const double *kx,
	const double *ky,
	std::complex<double> omega,
	const double &thickness, const double r[2],
	const std::complex<double> *q, // length 2*glist.n
	const std::complex<double> *kp, // size (2*glist.n)^2 (k-parallel matrix)
	const std::complex<double> *phi, // size (2*glist.n)^2
	const std::complex<double> *epsilon_inv, // size (glist.n)^2
	const std::complex<double> *epsilon2, // size (2*glist.n)^2
	int epstype,
	const std::complex<double> *ab, // length 4*glist.n
	double integral[6],
	std::complex<double> *work
){
	static const std::complex<double> zero(0.);
	static const std::complex<double> one(1.);
	const size_t n2 = 2*n;
	const size_t n4 = 2*n2;
	const std::complex<double> iomega(1./omega);

	std::complex<double> *f = work;
	if(NULL == work){
		f = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>)*12*n4);
	}
	std::complex<double> *g = f + 6*n4;
	RNP::TBLAS::Fill(6*n4, 0., f, 1);

	// We will simultaneously generate f for all 6 components

	// Generate the Fourier phase factors
	for(size_t i = 0; i < n; ++i){ // this is for ex
		double theta = (kx[i]*r[0] + ky[i]*r[1]);
		f[0*n4+1*n+i] = std::complex<double>(cos(theta),-sin(theta)); // note the minus
	}
	RNP::TBLAS::Copy(n, &f[n],1, &f[1*n4+0*n],1); // -ey
	if(EPSILON2_TYPE_BLKDIAG1_SCALAR == epstype || EPSILON2_TYPE_BLKDIAG2_SCALAR == epstype){
		RNP::TBLAS::Copy(n, &f[n4],1, &f[2*n4+n2],1); // ez
		RNP::TBLAS::Scale(n, iomega*epsilon_inv[0], &f[2*n4+n2],1);
	}else{
		RNP::TBLAS::MultMV<'C'>(n,n, iomega,epsilon_inv,n, &f[n4],1, zero,&f[2*n4+n2],1); // ez
	}
	RNP::TBLAS::Copy(n, &f[2*n4+n2],1, &f[2*n4+n2+n],1); // ez
	for(size_t i = 0; i < n; ++i){ // ez
		f[2*n4+n2+i] *= ky[i];
	}
	for(size_t i = 0; i < n; ++i){ // ez
		f[2*n4+n2+n+i] *= -kx[i];
	}
	RNP::TBLAS::Copy(n, &f[n4],1, &f[3*n4+n2],1); // hx
	RNP::TBLAS::Copy(n, &f[n4],1, &f[4*n4+n2+n],1); // hy
	for(size_t i = 0; i < n; ++i){ // -hz
		f[5*n4+0*n+i] = f[n4+i] * iomega * kx[i];
	}
	for(size_t i = 0; i < n; ++i){ // -hz
		f[5*n4+1*n+i] = f[n4+i] * iomega * ky[i];
	}

	// Form M^H C^H * f
	// M = [ kp.phi.inv(q)   -kp.phi.inv(q) ]
	//     [     phi                phi     ]
	if(NULL != phi){
		MultKPMatrix("C", omega, n, kx, ky, epsilon_inv, epstype, kp, 6, f,n4, &g[n2],n4);
		RNP::TBLAS::MultMM<'C','N'>(n2,6,n2, one,phi,n2, &g[n2],n4, zero,&g[0 ],n4);
		RNP::TBLAS::MultMM<'C','N'>(n2,6,n2, one,phi,n2, &f[n2],n4, zero,&g[n2],n4);
	}else{
		MultKPMatrix("C", omega, n, kx, ky, epsilon_inv, epstype, kp, 6, f,n4, &g[0],n4);
		RNP::TBLAS::CopyMatrix<'A'>(n2,6, &f[n2],n4, &g[n2],n4);
	} // top half of g is phi^H kp^H ftop, bottom half is phi^H fbot
	for(size_t i = 0; i < n2; ++i){
		RNP::TBLAS::Scale(6, 1./std::conj(omega*q[i]), &g[i],n4);
	}
	for(size_t j = 0; j < 6; ++j){
		for(size_t i = 0; i < n2; ++i){
			std::complex<double> t(g[i+j*n4]);
			g[i+j*n4] += g[i+n2+j*n4];
			g[i+n2+j*n4] -= t;
		}
	}

	for(size_t i = 0; i < 6; ++i){
		integral[i] = 0;
	}

	// At this point, each column of g contains a vector v such that
	// v v^H is a rank-1 matrix that must be integrated over.
	const double d_2 = 0.5*thickness;
	for(size_t j = 0; j < n4; ++j){
		const double sj = (j < n2) ? 1 : -1;
		const size_t qj = (j < n2) ? j : j-n2;
		for(size_t i = 0; i < n4; ++i){
			const double si = (i < n2) ? 1 : -1;
			const size_t qi = (i < n2) ? i : i-n2;
			std::complex<double> term = std::conj(ab[i]) * ab[j];
			term *= thickness*std::exp(std::complex<double>(0,d_2) * (q[qj] - std::conj(q[qi])));
			term *= zsinc(d_2 * (sj*q[qj] - si*std::conj(q[qi])));
			for(size_t ii = 0; ii < 6; ++ii){
				integral[ii] += std::real(std::conj(g[j+ii*n4])*g[i+ii*n4] * term);
			}
		}
	}

	if(NULL == work){
		rcwa_free(f);
	}
}
