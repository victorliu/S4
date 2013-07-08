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
# include <IO.h>
# include <cstdio>

#include <numalloc.h>

static inline void* rcwa_malloc(size_t size){
	void *ret = malloc_aligned(size, 16);
	//memset(ret, 0x0, size);
	return ret;
}
static inline void rcwa_free(void *ptr){
	free_aligned(ptr);
}

//typedef int integer;
extern "C" void zgelss_(
	const integer &m, const integer &n, const integer &nRHS,
	std::complex<double> *a, const integer &lda,
	std::complex<double> *b, const integer &ldb,
	double *s, const double &rcond, integer *rank,
	std::complex<double> *work, const integer &lwork,
	double *rwork, integer *info
);
static void SingularLinearSolve(
	size_t m, size_t n, size_t nRHS, std::complex<double> *a, size_t lda,
	std::complex<double> *b, size_t ldb, const double &rcond
){
	integer info, rank;
	std::complex<double> dummy;
	zgelss_(m, n, nRHS, a, lda, b, ldb, NULL, rcond, &rank, &dummy, -1, NULL, &info);
	integer lwork = (integer)dummy.real();
	std::complex<double> *work = (std::complex<double>*)rcwa_malloc(sizeof(std::complex<double>) * lwork);
	double *rwork = (double*)rcwa_malloc(sizeof(double) * 6*m);
	
	zgelss_(m, n, nRHS, a, lda, b, ldb, rwork, rcond, &rank, work, lwork, rwork+m, &info);
	
	rcwa_free(rwork);
	rcwa_free(work);
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

	if((size_t)-1 == lwork){
		RNP::Eigensystem(n2, NULL, n2, q, NULL, 1, phi, n2, work_, NULL, lwork);
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
		std::cerr << "Layer eigensystem returned info = " << info << std::endl;
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
		ab[i]    *= std::exp(iq*dz);
		ab[i+n2] *= std::exp(iq*(thickness-dz));
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
	
	for(size_t i = 0; i < n; ++i){
		double theta = (kx[i]*r[0] + ky[i]*r[1]);

		std::complex<double> phase(cos(theta),sin(theta));
		fH[0] += hx[i]*phase;
		fH[1] += hy[i]*phase;
		fE[0] += ex[i]*phase;
		fE[1] -= ney[i]*phase;
		fH[2] += (kx[i] * -ney[i] - ky[i] * ex[i]) / omega *phase;
		eh[i] = (ky[i]*hx[i] - kx[i]*hy[i]) *phase;
	}
	
	if(NULL != efield && NULL != epsilon_inv){
		efield[0] = fE[0];
		efield[1] = fE[1];
		efield[2] = 0;
		RNP::TBLAS::MultMV<'N'>(n,n, z_one,epsilon_inv,n, eh,1, z_zero,&eh[n],1);
		for(size_t i = 0; i < n; ++i){
			efield[2] += eh[n+i]/omega;
		}
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
		dz[i] = (hy[i]*hx[i] - kx[i]*hy[i]) / omega;
	}
	RNP::TBLAS::MultMV<'N'>(n,n, z_one,epsilon_inv,n, dz,1, z_zero,ez,1);
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
	double *integral,
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
		RNP::TBLAS::MultMM<'C','N'>(n,n,n, z_one,epsilon_inv,n, epsilon_inv,n, z_zero,&Q[n2+n2*n4],n4);
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
			*integral += term.real();
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
	RNP::TBLAS::MultMV<'C'>(n,n, iomega,epsilon_inv,n, &f[n4],1, zero,&f[2*n4+n2],1); // ez
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
